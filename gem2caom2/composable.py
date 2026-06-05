# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2026.                            (c) 2026.
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
import sys
import traceback

from caom2pipe.client_composable import ClientCollection
from caom2pipe import manage_composable as mc
from caom2pipe import run_composable as rc
from gem2caom2 import ghost_preview_augmentation, preview_augmentation
from gem2caom2 import pull_augmentation, data_source
from gem2caom2 import cleanup_augmentation, fits2caom2_augmentation
from gem2caom2 import svofps
from gem2caom2.gem_name import GemName
from gem2caom2.program_metadata import MDContext, PIMetadata


DATA_VISITORS = [ghost_preview_augmentation]
META_VISITORS = [pull_augmentation, fits2caom2_augmentation, preview_augmentation, cleanup_augmentation]


class GemClientCollection(ClientCollection):
    """
    Extend ClientCollection to have a place to hold and reference the
    archive.gemini.edu and svo sessions.
    """

    def __init__(self, config):
        super().__init__(config)
        self._gemini_session = None
        self._svo_session = None

    @property
    def gemini_session(self):
        return self._gemini_session

    @gemini_session.setter
    def gemini_session(self, value):
        self._gemini_session = value

    @property
    def svo_session(self):
        return self._svo_session

    @svo_session.setter
    def svo_session(self, value):
        self._svo_session = value


def _common_init():
    config = mc.Config()
    config.get_executors()
    clients = GemClientCollection(config)
    meta_visitors = META_VISITORS
    gemini_session = mc.get_endpoint_session()
    svofps_session = mc.get_endpoint_session()
    filter_cache = svofps.FilterMetadataCache(svofps_session)
    pi_metadata = PIMetadata(gemini_session)
    md_context = MDContext(filter_cache, pi_metadata)
    clients.gemini_session = gemini_session
    clients.svo_session = svofps_session
    if config.use_local_files or mc.TaskType.SCRAPE in config.task_types:
        meta_visitors = [
            fits2caom2_augmentation,
            preview_augmentation,
            cleanup_augmentation,
        ]
    mc.StorageName.collection = config.collection
    mc.StorageName.namespace = config.namespace
    mc.StorageName.scheme = config.scheme
    mc.StorageName.preview_scheme = config.preview_scheme
    return clients, config, meta_visitors, md_context


def _run():
    """
    Uses a todo file with file names, even though Gemini provides
    information about existing data referenced by observation ID.
    """
    clients, config, meta_visitors, md_context = _common_init()
    if config.use_local_files or mc.TaskType.SCRAPE in config.task_types:
        source = data_source.GeminiListDirSeparateDataSource(config, md_context)
    else:
        source = data_source.GeminiTodoFile(config, md_context)
    return rc.run_by_todo_runner_meta(
        config=config,
        meta_visitors=meta_visitors,
        data_visitors=DATA_VISITORS,
        sources=[source],
        clients=clients,
        organizer_module_name='gem2caom2.gemini_metadata',
        organizer_class_name='GeminiOrganizeExecutesRunnerMeta',
        storage_name_ctor=GemName,
    )


def run():
    """Wraps _run in exception handling, with sys.exit calls."""
    try:
        result = _run()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_by_public():
    """Run the processing for observations that are public, but there are
    no artifacts representing the previews in CAOM, or a FITS file in ad.

    Called as gem_run_public. The time-boxing is based on timestamps from a
    state.yml file. Call once/day, since data release timestamps have times
    of 00:00:00.000.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    clients, config, meta_visitors, md_context = _common_init()
    incremental_source = data_source.PublicIncremental(config, clients.query_client, md_context)
    return rc.run_by_state_runner_meta(
        config=config,
        meta_visitors=meta_visitors,
        data_visitors=DATA_VISITORS,
        sources=[incremental_source],
        clients=clients,
        organizer_module_name='gem2caom2.gemini_metadata',
        organizer_class_name='GeminiOrganizeExecutesRunnerMeta',
        storage_name_ctor=GemName,
    )


def run_by_public():
    try:
        result = _run_by_public()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_state():
    """Run incremental processing for observations that are posted on the site
    archive.gemini.edu. This depends on the incremental query endpoint.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    clients, config, meta_visitors, md_context = _common_init()
    incremental_source = data_source.IncrementalSource(config, clients.gemini_session, md_context)
    return rc.run_by_state_runner_meta(
        config=config,
        meta_visitors=meta_visitors,
        data_visitors=DATA_VISITORS,
        sources=[incremental_source],
        clients=clients,
        organizer_module_name='gem2caom2.gemini_metadata',
        organizer_class_name='GeminiOrganizeExecutesRunnerMeta',
        storage_name_ctor=GemName,
    )


def run_state():
    try:
        result = _run_state()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_incremental_diskfiles():
    """Run incremental processing for observations that are posted on the site archive.gemini.edu. This depends on
    the diskfiles query endpoint.

    :return 0 if successful, -1 if there's any sort of failure.
    """
    clients, config, meta_visitors, md_context = _common_init()
    incremental_source = data_source.IncrementalSourceDiskfiles(config, clients.gemini_session, GemName, md_context)
    return rc.run_by_state_runner_meta(
        config=config,
        meta_visitors=meta_visitors,
        data_visitors=DATA_VISITORS,
        sources=[incremental_source],
        clients=clients,
        organizer_module_name='gem2caom2.gemini_metadata',
        organizer_class_name='GeminiOrganizeExecutesRunnerMeta',
        storage_name_ctor=GemName,
    )


def run_incremental_diskfiles():
    try:
        result = _run_incremental_diskfiles()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
