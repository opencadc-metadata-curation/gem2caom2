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

from bs4 import BeautifulSoup
from collections import defaultdict, deque, OrderedDict
from datetime import datetime
from os import scandir

from cadcdata import FileInfo
from caom2utils.blueprints import _to_int
from caom2utils.data_util import get_file_type
from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe.manage_composable import CaomName, ISO_8601_FORMAT, make_datetime, query_endpoint_session
from gem2caom2.gem_name import GemName
from gem2caom2.obs_file_relationship import repair_data_label
from gem2caom2.util import set_instrument_case


__all__ = ['GEM_BOOKMARK', 'IncrementalSource', 'PublicIncremental']

GEM_BOOKMARK = 'gemini_timestamp'
MAX_ENTRIES = 500


class IncrementalSource(dsc.IncrementalDataSource):
    """Implements the identification of the work to be done, by querying
    archive.gemini.edu's incremental endpoint, in time-boxed chunks.

    OO - email 21-04-21
    entrytime is when the DB record behind the JSON being displayed was
    created.
    """

    def __init__(self, config, session, md_context):
        super().__init__(config, start_key=GEM_BOOKMARK)
        self._max_records_encountered = False
        self._session = session
        self.md_context = md_context

    def _initialize_end_dt(self):
        self._end_dt = datetime.now()

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        """
        :param prev_exec_dt datetime start of the time-boxed chunk
        :param exec_dt datetime end of the time-boxed chunk
        :return: a deque of StorageName instances with time their associated JSON (DB) records were modified from
            archive.gemini.edu.
        """
        self._logger.debug(f'Begin get_time_box_work from {prev_exec_dt} to {exec_dt}.')
        self._max_records_encountered = False
        # datetime format 2019-12-01T00:00:00.000000
        prev_dt_str = prev_exec_dt.strftime(ISO_8601_FORMAT)
        exec_dt_str = exec_dt.strftime(ISO_8601_FORMAT)
        url = (
            f'https://archive.gemini.edu/jsonsummary/canonical/'
            f'NotFail/notengineering/'
            f'entrytimedaterange={prev_dt_str}%20{exec_dt_str}/'
            f'?orderby=entrytime'
        )

        # needs to be ordered by timestamps when processed
        self._logger.info(f'Querying {url}')
        entries = deque()
        response = None
        try:
            response = query_endpoint_session(url, self._session)
            if response is None:
                self._logger.warning(f'Could not query {url}.')
            else:
                metadata = response.json()
                response.close()
                if metadata is not None:
                    if len(metadata) == 0:
                        self._logger.warning(
                            f'No query results returned for interval from {prev_exec_dt} to {exec_dt}.'
                        )
                    else:
                        for entry in metadata:
                            file_name = entry.get('name')
                            entry_dt = make_datetime(entry.get('entrytime'))
                            entries.append(
                                dsc.RunnerMeta(
                                    GemName(file_name=file_name, md_context=self.md_context), entry_dt
                                )
                            )
        finally:
            if response is not None:
                response.close()
        if len(entries) == MAX_ENTRIES:
            self._max_records_encountered = True
        self._reporter.capture_todo(len(entries), 0, 0)
        self._logger.debug('End get_time_box_work.')
        return entries

    def max_records_encountered(self):
        return self._max_records_encountered


class PublicIncremental(dsc.QueryTimeBoxDataSource):
    """Implements the identification of the work to be done, by querying
    the local TAP service for files that have recently gone public."""

    def __init__(self, config, query_client, md_context):
        super().__init__(config)
        self._query_client = query_client
        self.md_context = md_context

    def _initialize_end_dt(self):
        self._end_dt = datetime.now()

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        """
        :param prev_exec_dt datetime start of the timestamp chunk
        :param exec_dt datetime end of the timestamp chunk
        :return: a list of file names with time they were modified in /ams,
            structured as an astropy Table (for now).
        """
        self._logger.debug('Begin get_time_box_work')
        # datetime format 2019-12-01T00:00:00.000000
        prev_dt_str = prev_exec_dt.strftime(ISO_8601_FORMAT)
        exec_dt_str = exec_dt.strftime(ISO_8601_FORMAT)
        query = f"""SELECT O.observationID, A.uri, A.lastModified
            FROM caom2.Observation AS O
            JOIN caom2.Plane ON O.obsID = caom2.Plane.obsID
            JOIN caom2.Artifact AS A ON caom2.Plane.planeID = A.planeID
            WHERE O.collection = '{self._config.collection}'
            AND caom2.Plane.dataRelease > '{prev_dt_str}'
            AND caom2.Plane.dataRelease <= '{exec_dt_str}'
            AND (
                SELECT COUNT(*)
                FROM caom2.Artifact AS A2
                WHERE A2.planeID = caom2.Plane.planeID
            ) = 1
            ORDER BY O.maxLastModified ASC
        """
        result = clc.query_tap_client(query, self._query_client)
        # results look like:
        # GN-2019B-ENG-1-160-008, gemini:GEM/N20191202S0125.fits, ISO 8601

        entries = deque()
        for row in result:
            gem_name = GemName(file_name=CaomName(row['uri']).file_name, md_context=self.md_context)
            gem_name._obs_id = row['observationID']
            entries.append(dsc.RunnerMeta(gem_name, make_datetime(row['lastModified'])))
        self._reporter.capture_todo(len(entries), 0, 0)
        self._logger.debug('End get_time_box_work')
        return entries


class IncrementalSourceDiskfiles(dsc.IncrementalDataSource):
    """Implements the identification of the work to be done, by querying archive.gemini.edu's incremental diskfiles
    endpoint, in time-boxed chunks.

    OO - email 21-04-21
    entrytime is when the DB record behind the JSON being displayed was created.
    """

    def __init__(self, config, gemini_session, storage_name_ctor, md_context):
        super().__init__(config, start_key=GEM_BOOKMARK)
        self._max_records_encountered = False
        self._session = gemini_session
        self._storage_name_ctor = storage_name_ctor
        self.md_context = md_context

    def _initialize_end_dt(self):
        self._end_dt = datetime.now()

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        """
        :param prev_exec_dt datetime start of the time-boxed chunk
        :param exec_dt datetime end of the time-boxed chunk
        :return: a deque of file names with time their associated JSON (DB) records were modified from
            archive.gemini.edu.
        """

        self._logger.debug(f'Begin get_time_box_work from {prev_exec_dt} to {exec_dt}.')
        self._max_records_encountered = False
        # datetime format 2019-12-01T00:00:00 => no microseconds in the url
        prev_exec_dt_iso = prev_exec_dt.replace(microsecond=0)
        exec_dt_iso = exec_dt.replace(microsecond=0)
        url = (
            f'https://archive.gemini.edu/diskfiles/NotFail/notengineering/not_site_monitoring/canonical/entrytimedaterange='
            f'{prev_exec_dt_iso.isoformat()}--{exec_dt_iso.isoformat()}'
        )

        # needs to be ordered by timestamps when processed
        self._logger.info(f'Querying {url}')
        entries = deque()
        response = None
        try:
            response = query_endpoint_session(url, self._session)
            if response is None:
                self._logger.warning(f'Could not query {url}.')
            else:
                # x is a dict, sorted by the key
                # key is the timestamp
                # value is the html meta that comes back, also in a dict
                metadata = self._parse_diskfiles_response(response.text)
                response.close()
                if len(metadata) == 0:
                    self._logger.warning(f'No query results returned for interval from {prev_exec_dt} to {exec_dt}.')
                else:
                    for entry_dt, values in metadata.items():
                        for value in values:
                            file_name = value.get('filename')
                            storage_name = self._storage_name_ctor(file_name, self.md_context)
                            storage_name.file_info[storage_name.destination_uris[0]] = FileInfo(
                                id=file_name,
                                size=_to_int(value.get('data_size')),
                                md5sum=value.get('data_md5'),
                                file_type=get_file_type(file_name),
                            )
                            repaired_data_label = repair_data_label(
                                file_name, value.get('datalabel'), set_instrument_case(value.get('instrument'))
                            )
                            storage_name.obs_id = repaired_data_label
                            storage_name._fullheader = value.get('fullheader')
                            entries.append(dsc.RunnerMeta(storage_name, entry_dt))
        finally:
            if response is not None:
                response.close()
        if len(entries) == MAX_ENTRIES:
            self._max_records_encountered = True
        self._reporter.capture_todo(len(entries), 0, 0)
        self._logger.debug('End get_time_box_work.')
        return entries

    def _parse_diskfiles_response(self, html_string):
        temp = defaultdict(list)
        soup = BeautifulSoup(html_string, features='lxml')
        rows = soup.find_all('tr', {'class': 'alternating'})
        for row in rows:
            cells = row.find_all('td')
            entry_time = make_datetime(cells[5].text.strip())  # what the https query is keyed on
            file_name = cells[0].text.strip()
            if file_name.endswith('-fits!-md!'):
                continue
            fullheader = cells[0].find_all('a')[0].get('href')
            instrument = cells[3].find_all('a')[0].get('href')
            value = {
                'filename': file_name,
                'datalabel': cells[1].text.strip().split('\n')[-1],
                'instrument': instrument.split('/')[-1],
                'lastmod': make_datetime(cells[6].text.strip()),
                'data_size': cells[10].text.strip(),
                'data_md5': cells[11].text.strip(),
                'fullheader': fullheader.split('/')[-1],
            }
            temp[entry_time].append(value)
        result = OrderedDict(sorted(temp.items()))
        return result

    def max_records_encountered(self):
        return self._max_records_encountered


class GeminiTodoFile(dsc.TodoFileDataSourceRunnerMeta):

    def __init__(self, config, md_context):
        super().__init__(config, GemName)
        self.md_context = md_context

    def _find_work(self, entry_path):
        with open(entry_path) as f:
            for line in f:
                temp = line.strip()
                if len(temp) > 0:
                    # ignore empty lines
                    self._logger.debug(f'Adding entry {temp} to work list.')
                    self._work.append(GemName(file_name=temp, md_context=self.md_context))


class GeminiListDirSeparateDataSource(dsc.ListDirSeparateDataSource):

    def __init__(self, config, md_context):
        super().__init__(config)
        self.md_context = md_context

    def _append_work(self, entry):
        with scandir(entry) as dir_listing:
            for entry in dir_listing:
                if entry.is_dir() and self._recursive:
                    self._append_work(entry.path)
                else:
                    if self.default_filter(entry):
                        self._logger.debug(f'Adding {entry.path} to work list.')
                        self._work.append(GemName(source_names=[entry.path], md_context=self.md_context))
