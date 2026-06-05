# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

import logging
import traceback

from datetime import datetime, timezone
from os.path import basename, exists, join

import matplotlib.image as image

from cadcutils import exceptions
from caom2 import Observation, ProductType, ReleaseType
from caom2pipe.caom_composable import update_file_info
from caom2pipe import manage_composable as mc
from gem2caom2.util import Inst

__all__ = ['visit']


PREVIEW_URL = 'https://archive.gemini.edu/preview/'
MIME_TYPE = 'image/jpeg'


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    working_dir = kwargs.get('working_directory', './')
    clients = kwargs.get('clients')
    if clients is None or clients.data_client is None:
        logging.warning('Need a cadc_client to retrieve preview records.')
        return observation
    reporter = kwargs.get('reporter')
    if reporter is None:
        raise mc.CadcException('Visitor needs a reporter parameter.')
    observable = reporter._observable
    storage_name = kwargs.get('storage_name')
    if storage_name is None:
        raise mc.CadcException('Visitor needs a storage_name parameter.')

    if observation.instrument.name == Inst.GHOST.value:
        logging.info(f'Skip generic preview augmentation for GHOST.')
        return observation

    count = 0
    for plane in observation.planes.values():
        if (
            plane.data_release is None
            # data_release is timezone naive
            or plane.data_release > datetime.now(tz=timezone.utc).replace(tzinfo=None)
        ):
            logging.info(f'Plane {plane.product_id} is proprietary. No preview access or thumbnail creation.')
            continue
        if plane.product_id != storage_name.product_id:
            continue
        count += _do_prev(
            observation.observation_id,
            working_dir,
            plane,
            clients,
            observable,
            storage_name,
        )
    result = {'files': count}
    logging.info(f'Completed preview augmentation for {observation.observation_id}. {count} files modified.')
    return observation


def _do_prev(obs_id, working_dir, plane, clients, observable, gem_name):
    """Retrieve the preview file, so that a thumbnail can be made,
    store the preview if necessary, and the thumbnail, to CADC storage.
    Then augment the CAOM observation with the two additional artifacts.
    """
    logging.debug('Begin _do_prev')
    count = 0

    if observable.rejected.is_no_preview(gem_name.prev):
        logging.info(f'Stopping visit because no preview exists for {gem_name.prev} in observation {obs_id}.')
        observable.rejected.record(mc.Rejected.NO_PREVIEW, gem_name.prev)
        count += _check_for_delete(gem_name.prev, gem_name.prev_uri, observable, plane)
    else:
        count = get_representation(clients, gem_name, working_dir, plane, observable)
    logging.debug('End _do_prev')
    return count


def _check_for_delete(file_name, uri, observable, plane):
    """If the preview file doesn't exist, but the artifact that represents it
    does, remove that artifact from the Observation instance."""
    result = 0
    if observable.rejected.is_no_preview(file_name) and uri in plane.artifacts.keys():
        logging.warning(f'Removing artifact for non-existent preview {uri}')
        plane.artifacts.pop(uri)
        result = 1
    return result


def _augment(plane, uri, fqn, product_type):
    count = 0
    temp = None
    if uri in plane.artifacts:
        temp = plane.artifacts[uri]
    plane.artifacts[uri] = mc.get_artifact_metadata(fqn, product_type, ReleaseType.DATA, uri, temp)
    return count


def _retrieve_from_gemini(
    gem_name,
    observable,
    plane,
    preview_fqn,
):
    temp = basename(gem_name.file_name)
    preview_url = f'{PREVIEW_URL}{temp}'
    logging.debug(f'Retrieve {gem_name.prev} from archive.gemini.edu at {preview_url}.')
    try:
        mc.http_get(preview_url, preview_fqn)
    except Exception as e:
        logging.warning(e)
        if observable.rejected.check_and_record(str(e), gem_name.prev):
            _check_for_delete(gem_name.prev, gem_name.prev_uri, observable, plane)


def get_representation(clients, storage_name, working_dir, plane, observable):
    logging.debug(f'Begin get_representation for {storage_name.file_uri}')
    count = 0
    prev_info = clients.data_client.info(storage_name.prev_uri)
    thumb_info = clients.data_client.info(storage_name.thumb_uri)

    if prev_info and thumb_info:
        logging.info(f'Updating FileInfo for {storage_name.file_uri}')
        update_file_info(prev_info, storage_name.prev_uri, plane)
        update_file_info(thumb_info, storage_name.thumb_uri, plane)
    else:
        preview_fqn = join(working_dir, storage_name.prev)
        thumb_fqn = join(working_dir, storage_name.thumb)
        if prev_info:
            # get the preview from CADC
            logging.debug(f'Retrieve {storage_name.prev} from CADC at {storage_name.prev_uri}.')
            clients.data_client.get(working_dir, storage_name.prev_uri)
        else:
            # get the preview from Gemini
            _retrieve_from_gemini(storage_name, observable, plane, preview_fqn)
            if exists(preview_fqn):
                clients.data_client.put(working_dir, storage_name.prev_uri)
                count += 1
            else:
                logging.error(f'Could not retrieve {storage_name.prev_uri} from archive.gemini.edu')
                return count

        # build the thumbnail - include the retry for failure to pull down the preview from Gemini
        try:
            image.thumbnail(preview_fqn, thumb_fqn, scale=0.25)
        except ValueError as _:
            # probably the jpg did not transfer properly, so try to retrieve it one more time
            #
            # have a retry here, because otherwise there's no way to update the file in CADC storage without
            # intervention from Ops - i.e. the file will retrieve from CADC, so there will be no succeeding
            # attempt to retrieve from Gemini that might otherwise fix the value
            logging.debug(traceback.format_exc())
            logging.warning(
                f'matplotlib error handling {storage_name.prev}.Try to retrieve from {PREVIEW_URL} one more time.'
            )
            _retrieve_from_gemini(storage_name, observable, plane, preview_fqn)
            if exists(preview_fqn):
                clients.data_client.put(working_dir, storage_name.prev_uri)
                count += 1
                image.thumbnail(preview_fqn, thumb_fqn, scale=0.25)
            else:
                return count

        clients.data_client.put(working_dir, storage_name.thumb_uri)
        count += 1

        _augment(plane, storage_name.prev_uri, preview_fqn, ProductType.PREVIEW)
        _augment(plane, storage_name.thumb_uri, thumb_fqn, ProductType.THUMBNAIL)
    logging.debug('End get_representation')
    return count
