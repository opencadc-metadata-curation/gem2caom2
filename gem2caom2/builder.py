# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from caom2pipe import name_builder_composable as nbc
from gem2caom2 import gem_name, external_metadata


__all__ = ['EduQueryBuilder', 'GemObsIDBuilder', 'get_instrument']


class EduQueryBuilder(nbc.StorageNameBuilder):
    """
    Get the file metadata by querying archive.gemini.edu. This information is
    required to find the data label for a file name, so that a StorageName
    instance can be built.

    This class delays the time when the metadata must be queried until
    just before it is used by the execute_composable methods, so that the
    memory footprint of the pipeline does not have to support the
    gemini-sourced metadata of all entries in the list of work to be done.
    """

    def __init__(self, config):
        super(EduQueryBuilder, self).__init__()
        self._todo_list = None

    @property
    def todo_list(self):
        return self._todo_list

    @todo_list.setter
    def todo_list(self, to_list):
        self._todo_list = {value: key for key, value in to_list.items()}

    def build(self, entry):
        """
        :param entry: a Gemini file name
        :return: an instance of StorageName for use in execute_composable.
        """
        if self._config.use_local_files:
            raise NotImplementedError('The need has not been encountered '
                                      'in the real world.')

        external_metadata.get_obs_metadata(
            gem_name.GemName.remove_extensions(entry))
        instrument = get_instrument()
        storage_name = gem_name.GemName(file_name=entry, instrument=instrument)
        return storage_name


class GemObsIDBuilder(nbc.StorageNameBuilder):
    """
    To be able to build a StorageName instance with an observation ID.
    """

    def __init__(self, config):
        super(GemObsIDBuilder, self).__init__()
        self._config = config
        self._instrument = None
        self._logger = logging.getLogger(__name__)

    def _read_instrument_locally(self, entry):
        self._logger.debug(f'Use a local file to read instrument from the '
                           f'headers for {entry}.')
        headers = ac.read_fits_headers(
            f'{self._config.working_directory}/{entry}')
        self._instrument = external_metadata.Inst(headers[0].get('INSTRUME'))

    def _read_instrument_remotely(self, entry):
        self._logger.debug(
            'Read instrument from archive.gemini.edu.')
        file_id = gem_name.GemName.remove_extensions(entry)
        external_metadata.get_obs_metadata(file_id)
        self._instrument = get_instrument()

    def build(self, entry):
        """
        :param entry: a Gemini file name or observation ID, depending on
            the configuration
        :return: an instance of StorageName for use in execute_composable.
        """
        self._logger.debug(f'Build a StorageName instance for {entry}.')
        try:
            if self._config.features.supports_latest_client:
                if (mc.TaskType.SCRAPE in self._config.task_types or
                        self._config.use_local_files):
                    self._read_instrument_locally(entry)
                    result = gem_name.GemName(file_name=entry,
                                              instrument=self._instrument,
                                              v_collection=gem_name.COLLECTION,
                                              v_scheme=gem_name.V_SCHEME,
                                              entry=entry)
                elif self._config.features.use_file_names:
                    self._read_instrument_remotely(entry)
                    result = gem_name.GemName(file_name=entry,
                                              instrument=self._instrument,
                                              v_collection=gem_name.COLLECTION,
                                              v_scheme=gem_name.V_SCHEME,
                                              entry=entry)
                else:
                    raise mc.CadcException('The need has not been encountered '
                                           'in the real world yet.')
            else:
                if (mc.TaskType.INGEST_OBS in self._config.task_types and
                        '.fits' not in entry):
                    # anything that is NOT ALOPEKE/ZORRO, which are the only
                    # two instruments that change the behaviour of the
                    # GemName constructor - and yeah, that abstraction is
                    # leaking like a sieve.
                    self._logger.debug('INGEST_OBS, hard-coded instrument.')
                    instrument = external_metadata.Inst.CIRPASS
                    result = gem_name.GemName(obs_id=entry,
                                              instrument=instrument,
                                              entry=entry)
                elif (mc.TaskType.SCRAPE in self._config.task_types or
                        self._config.use_local_files):
                    self._read_instrument_locally(entry)
                    result = gem_name.GemName(file_name=entry,
                                              instrument=self._instrument,
                                              entry=entry)
                elif self._config.features.use_file_names:
                    self._read_instrument_remotely(entry)
                    result = gem_name.GemName(file_name=entry,
                                              instrument=self._instrument,
                                              entry=entry)
                else:
                    raise mc.CadcException('The need has not been encountered '
                                           'in the real world yet.')
            self._logger.debug('Done build.')
            return result
        except Exception as e:
            self._logger.error(e)
            self._logger.debug(traceback.format_exc())
            raise mc.CadcException(e)


def get_instrument():
    inst = external_metadata.om.get('instrument')
    if inst == 'ALOPEKE':
        # because the value in JSON is a different case than the value in
        # the FITS header
        inst = 'Alopeke'
    if inst == 'ZORRO':
        inst = 'Zorro'
    return external_metadata.Inst(inst)
