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

import json
import logging
import os
import shutil
import traceback

from astropy.io.fits import Header
from collections import deque
from datetime import datetime, timedelta, timezone
from traceback import format_exc
from unittest.mock import ANY, call, patch, Mock, PropertyMock
import gem_mocks

from cadcdata import FileInfo
from caom2 import SimpleObservation, Algorithm
from caom2pipe.data_source_composable import RunnerMeta
from caom2pipe.manage_composable import Config, make_datetime, State, TaskType, write_as_yaml
from gem2caom2 import composable, gem_name
from gem2caom2.data_source import GEM_BOOKMARK


STATE_FILE = f'{gem_mocks.TEST_DATA_DIR}/state.yml'
TODO_FILE = f'{gem_mocks.TEST_DATA_DIR}/todo.txt'
REJECTED_FILE = f'{gem_mocks.TEST_DATA_DIR}/logs/rejected.yml'
PROGRESS_FILE = f'{gem_mocks.TEST_DATA_DIR}/logs/progress.txt'
PUBLIC_TEST_JSON = f'{gem_mocks.TEST_DATA_DIR}/json/GN-2019B-ENG-1-160-008.json'


@patch('gem2caom2.composable.GemClientCollection')
@patch('caom2pipe.execute_composable.OrganizeExecutesRunnerMeta.do_one')
def test_run(run_mock, clients_mock, test_config, tmp_path, change_test_dir):
    # use a todo.txt file to drive work
    test_f_id = 'S20070130S0048'
    test_f_name = f'{test_f_id}.fits'
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.INGEST]
    gem_mocks.use_test_data_cache(test_config)
    test_config.write_to_file(test_config)

    with open(test_config.work_fqn, 'w') as f:
        f.write(f'{test_f_name}\n')

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    run_mock.return_value = (0, None)

    # execution
    composable._run()
    assert run_mock.called, 'should have been called'
    args, _ = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, gem_name.GemName), type(test_storage)
    # don't check obs_id, because it will be set by the do_one call in this scenario
    assert test_storage.file_name == test_f_name, 'wrong file name'
    clients_mock.return_value.metadata_client.read.assert_not_called(), 'meta read'


@patch('gem2caom2.composable.GemClientCollection')
@patch('caom2pipe.execute_composable.OrganizeExecutesRunnerMeta.do_one')
@patch('gem2caom2.data_source.query_endpoint_session')
@patch('caom2pipe.client_composable.query_tap_client')
def test_run_incremental_rc(
    tap_mock,
    query_mock,
    run_mock,
    clients_mock,
    test_config,
    tmp_path,
    change_test_dir,
):
    # use the original incremental endpoint to drive work
    query_mock.side_effect = gem_mocks.mock_query_endpoint_2
    tap_mock.side_effect = gem_mocks.mock_query_tap

    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'testproxy.pem'
    test_config.task_types = [TaskType.INGEST]
    test_config.interval = 600
    gem_mocks.use_test_data_cache(test_config)
    Config.write_to_file(test_config)
    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    _write_state(
        prior_timestamp='2021-01-01 20:03:00.000000',
        end_timestamp=datetime(year=2021, month=1, day=4, hour=23, minute=3, second=0),
        fqn=test_config.state_fqn,
    )
    run_mock.return_value = (0, None)
    composable._run_state()
    assert run_mock.called, 'run_mock should have been called'
    args, _ = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, gem_name.GemName), type(test_storage)
    # don't check obs_id, because it will be set by the do_one call in this scenario
    test_fid = 'N20210101S0042'
    assert test_storage.file_name == f'{test_fid}.fits', 'wrong file_name'
    assert test_storage.file_id == f'{test_fid}', 'wrong file_id'


@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
@patch('gem2caom2.composable.GemClientCollection')
@patch('caom2pipe.manage_composable.read_obs_from_file')
@patch('gem2caom2.data_source.query_endpoint_session')
def test_run_by_incremental2(
    query_mock,
    read_mock,
    clients_mock,
    exec_mock,
):
    clients_mock.data_client.return_value.info.side_effect = gem_mocks.mock_get_file_info
    clients_mock.data_client.return_value.get.side_effect = Mock()
    exec_mock.return_value = (0, None)
    clients_mock.metadata_client.create.side_effect = gem_mocks.mock_repo_create
    clients_mock.metadata_client.side_effect = gem_mocks.mock_repo_read
    clients_mock.metadata_client.update.side_effect = gem_mocks.mock_repo_update

    def _read_mock(ignore_fqn):
        return SimpleObservation(collection='TEST', observation_id='TEST_OBS_ID', algorithm=Algorithm('exposure'))

    read_mock.side_effect = _read_mock

    def _query_mock():
        x = json.loads(
            """[{
        "name": "2002feb11_0180.fits",
        "filename": "2002feb11_0180.fits.bz2",
        "path": "",
        "compressed": true,
        "file_size": 397169,
        "data_size": 1056960,
        "file_md5": "850d2df8b71c894984b4ff9d89cebcd8",
        "data_md5": "cf3ccf17539b23c4aa68c80308f885e2",
        "lastmod": "2019-12-17 00:21:08.159934+00:00",
        "mdready": true,
        "entrytime": "2019-12-17 00:21:08.178516+00:00",
        "size": 397169,
        "md5": "850d2df8b71c894984b4ff9d89cebcd8",
        "program_id": "GS-2002A-Q-8",
        "engineering": false,
        "science_verification": false,
        "calibration_program": false,
        "observation_id": "GS-2002A-Q-8-4",
        "data_label": "GS-2002A-Q-8-4-0180",
        "telescope": "Gemini-South",
        "instrument": "PHOENIX",
        "ut_datetime": "2002-04-04 21:47:53",
        "local_time": null,
        "observation_type": null,
        "observation_class": null,
        "object": "SZ Cha",
        "ra": 164.574204166667,
        "dec": -77.287963888889,
        "azimuth": null,
        "elevation": null,
        "cass_rotator_pa": null,
        "airmass": 1.478,
        "filter_name": "2150_(2)",
        "exposure_time": 240.0,
        "disperser": null,
        "camera": null,
        "central_wavelength": null,
        "wavelength_band": null,
        "focal_plane_mask": "107u_1.0-5.0 (8)",
        "detector_binning": "1x1",
        "detector_gain_setting": null,
        "detector_roi_setting": "Fixed",
        "detector_readspeed_setting": null,
        "detector_welldepth_setting": null,
        "detector_readmode_setting": "None",
        "spectroscopy": true,
        "mode": "spectroscopy",
        "adaptive_optics": false,
        "laser_guide_star": false,
        "wavefront_sensor": null,
        "gcal_lamp": null,
        "raw_iq": 50,
        "raw_cc": 50,
        "raw_wv": 80,
        "raw_bg": 80,
        "requested_iq": null,
        "requested_cc": null,
        "requested_wv": null,
        "requested_bg": null,
        "qa_state": "Undefined",
        "release": "2003-08-11",
        "reduction": "RAW",
        "types": "{'UNPREPARED', 'SOUTH', 'RAW', 'PHOENIX', 'GEMINI', 'SPECT'}",
        "phot_standard": null,
        "results_truncated": true
    }]"""
        )
        return x

    query_result = gem_mocks.Object()
    query_result.json = _query_mock
    query_mock.return_value = query_result

    exec_mock.return_value = (0, None)

    _write_cert()
    prior_s = datetime.now(tz=timezone.utc).timestamp() - 60
    _write_state(prior_s, fqn=f'{gem_mocks.TEST_DATA_DIR}/edu_query/state.yml')
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=f'{gem_mocks.TEST_DATA_DIR}/edu_query')
    try:
        # execution
        test_result = composable._run_state()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig

    assert exec_mock.called, 'exec mock not called'
    assert not (
        clients_mock.data_client.return_value.info.called
    ), 'data client mock get file info should not be not called'
    assert query_mock.called, 'query mock not called'


@patch('caom2pipe.client_composable.query_tap_client')
@patch('gem2caom2.composable.GemClientCollection')
@patch('gem2caom2.pull_augmentation.retrieve_headers')
@patch('gem2caom2.gemini_metadata.retrieve_json')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run_by_public(exec_mock, json_mock, header_mock, clients_mock, query_mock):
    exec_mock.side_effect = Mock(return_value=0)
    query_mock.side_effect = gem_mocks.mock_query_tap
    json_mock.side_effect = gem_mocks.mock_retrieve_json
    header_mock.side_effect = gem_mocks._mock_retrieve_headers
    expected_fqn = f'{gem_mocks.TEST_DATA_DIR}/logs/{gem_mocks.TEST_BUILDER_OBS_ID}.expected.xml'
    if not os.path.exists(expected_fqn):
        shutil.copy(f'{gem_mocks.TEST_DATA_DIR}/expected.xml', expected_fqn)

    _write_cert()
    now_dt = datetime.now(tz=timezone.utc).replace(tzinfo=None)
    prior_s = now_dt.timestamp() - 1440 * 60
    _write_state(prior_s)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=f'{gem_mocks.TEST_DATA_DIR}/edu_query')
    test_f_id = 'N20191101S0007'
    try:
        with patch(
            'caom2pipe.data_source_composable.QueryTimeBoxDataSource.end_dt', PropertyMock(return_value=now_dt)
        ):
            # execution
            test_result = composable._run_by_public()
            assert test_result == 0, 'wrong result'
    except Exception as e:
        logging.error(traceback.format_exc())
    finally:
        os.getcwd = getcwd_orig

    assert query_mock.called, 'tap mock not called'
    assert exec_mock.called, 'exec mock not called'
    args, _ = exec_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, gem_name.GemName), type(test_storage)
    assert test_storage.obs_id == 'GN-2019B-ENG-1-160-008', 'wrong obs id'
    assert test_storage.file_name == f'{test_f_id}.fits', 'wrong file_name'
    assert test_storage.file_id == test_f_id, 'wrong file_id'


@patch('gem2caom2.gemini_metadata.retrieve_json')
@patch('gem2caom2.pull_augmentation.http_get')
@patch('gem2caom2.svofps.FilterMetadataCache.filter_metadata')
@patch('gem2caom2.program_metadata.PIMetadata.get')
@patch('caom2pipe.client_composable.ClientCollection.metadata_client')
@patch('caom2pipe.client_composable.ClientCollection.data_client')
@patch('gem2caom2.pull_augmentation.retrieve_headers')
@patch('gem2caom2.data_source.query_endpoint_session')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
def test_run_by_incremental_reproduce(
    access_mock,
    query_mock,
    header_mock,
    data_client_mock,
    meta_client_mock,
    pi_mock,
    svo_mock,
    http_get_mock,
    json_mock,
    test_config,
    tmp_path,
    change_test_dir,
):
    # https://archive.gemini.edu/jsonsummary/canonical/NotFail/notengineering/
    # entrytimedaterange=
    # 2022-03-14T17:30:05.000006%202022-03-14T17:31:05.000006/
    # ?orderby=entrytime
    # get results
    query_mock.side_effect = gem_mocks.mock_query_endpoint_reproduce
    access_mock.return_value = 'https://localhost:2022'

    test_header = Header()
    test_header['INSTRUME'] = 'GMOS-S'
    header_mock.return_value = [test_header]
    data_client_mock.get_head.return_value = [test_header]
    meta_client_mock.read.return_value = None
    pi_mock.return_value = {}
    svo_mock.return_value = None
    json_mock.side_effect = gem_mocks.mock_retrieve_json

    def _repo_create_mock(observation):
        plane_count = 0
        artifact_count = 0
        for plane in observation.planes.values():
            plane_count += 1
            for _ in plane.artifacts.values():
                artifact_count += 1

        assert plane_count == 1, 'wrong plane count'
        assert artifact_count == 1, 'wrong artifact count'

    meta_client_mock.create = _repo_create_mock

    test_config.change_working_directory(tmp_path)
    test_config.logging_level = 'INFO'
    test_config.proxy_file_name = 'cadcproxy.pem'
    test_config.task_types = [TaskType.VISIT]
    test_config.interval = 70
    gem_mocks.use_test_data_cache(test_config)
    test_config.write_to_file(test_config)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    test_bookmark = {
        'bookmarks': {
            GEM_BOOKMARK: {
                'last_record': datetime.now() - timedelta(hours=1),
            },
        },
    }
    write_as_yaml(test_bookmark, test_config.state_fqn)

    # execution
    composable._run_state()
    assert meta_client_mock.read.called, 'should have been called'
    assert meta_client_mock.read.call_count == 2, f'wrong call count {meta_client_mock.read.call_count}'
    meta_client_mock.read.assert_called_with('GEMINI', 'test_data_label'), 'wrong run args'


@patch('gem2caom2.pull_augmentation.retrieve_headers')
@patch('gem2caom2.gemini_metadata.retrieve_json')
@patch('gem2caom2.composable.GemClientCollection')
@patch('gem2caom2.data_source.IncrementalSource.get_time_box_work', autospec=True)
def test_run_state_compression_commands(
    get_work_mock,
    clients_mock,
    json_mock,
    headers_mock,
    test_config,
    tmp_path,
    change_test_dir,
):
    test_config.change_working_directory(tmp_path)

    # this test works with FITS files, not header-only versions of FITS files, because it's testing the
    # decompression/recompression cycle but it's checking that the commands to the exec_cmd_array call are correct

    json_mock.side_effect = gem_mocks.mock_retrieve_json
    headers_mock.side_effect = gem_mocks._mock_retrieve_headers

    uris = {
        'GS-2005B-SV-301-16-005': FileInfo(
            'gemini:GEMINI/S20050825S0143.fits',
            size=19186560,  # not the compressed size of 4795130
            file_type='application/fits',
            md5sum='md5:24cf5c193a312d9aa76d94a5e2cf39c3',
        ),
    }

    def _mock_dir_list(arg1, output_file='', data_only=True, response_format='arg4'):
        result = deque()
        test_gem_name = gem_name.GemName(file_name='/test_files/S20050825S0143.fits.bz2', md_context=Mock())
        result.append(RunnerMeta(test_gem_name, datetime(2019, 10, 23, 16, 19)))
        return result

    get_work_mock.side_effect = _mock_dir_list
    clients_mock.return_value.data_client.info.side_effect = uris.get('GS-2005B-SV-301-16-005')

    test_config.task_types = [TaskType.STORE]
    test_config.logging_level = 'DEBUG'
    test_config.proxy_file_name = 'cadcproxy.pem'
    test_config.features.supports_latest_client = True
    test_config.features.supports_decompression = True
    test_config.use_local_files = True
    test_config.data_sources = ['/test_files']
    test_config.retry_failures = False
    test_config.cleanup_files_when_storing = False
    gem_mocks.use_test_data_cache(test_config)
    test_config.write_to_file(test_config)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    start_time = datetime.now() - timedelta(minutes=5)
    start_file_content = f'bookmarks:\n  gemini_timestamp:\n    last_record: {start_time}\n'
    with open(test_config.state_fqn, 'w') as f:
        f.write(start_file_content)

    try:
        test_result = composable._run_state()
        assert test_result == 0, 'expecting correct execution'
    except Exception as e:
        logging.error(e)
        logging.error(format_exc())
        raise e

    clients_mock.return_value.data_client.put.assert_called(), 'put'
    assert (
        clients_mock.return_value.data_client.put.call_count == 1
    ), 'put call count, no previews because it is just a STORE task'
    put_calls = [
        call(
            f'{tmp_path.as_posix()}/S20050825S0143',
            f'{test_config.scheme}:{test_config.collection}/S20050825S0143.fits',
        ),
    ]
    clients_mock.return_value.data_client.put.assert_has_calls(put_calls), 'wrong put args'

    # LocalStore, put is mocked, no info calls as part of that
    clients_mock.return_value.data_client.info.assert_not_called(), 'info'

    # LocalStore, get_head should not be called
    clients_mock.return_value.data_client.get_head.assert_not_called()
    # LocalStore, get should not be called
    clients_mock.return_value.data_client.get.assert_not_called()
    assert not clients_mock.return_value.metadata_client.read.called, 'read'


@patch('caom2pipe.manage_composable.ExecutionSummary', autospec=True)
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
def test_run_is_valid_fails(cap_mock, summary_mock, test_config, tmp_path):
    summary_mock.return_value.report.return_value = 'msg'
    cap_mock.return_value = 'https://localhost'
    test_f_id = ' N20220601S0052_ql_image'
    test_f_name = f'{test_f_id}.fits'
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.INGEST]
    orig_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        gem_mocks.use_test_data_cache(test_config)
        test_config.write_to_file(test_config)

        with open(test_config.work_fqn, 'w') as f:
            f.write(f'{test_f_name}\n')

        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')

        # execution
        test_result = composable._run()
        assert test_result == -1, 'expect failure'
        assert summary_mock.return_value.add_errors.called
        assert summary_mock.return_value.add_errors.call_count == 1
        assert summary_mock.return_value.add_rejections.called
        assert summary_mock.return_value.add_rejections.call_count == 2, 'one for capture_todo, one for rejection'
        assert not summary_mock.return_value.add_successes.called
        assert summary_mock.return_value.add_entries.called
        summary_mock.return_value.add_entries.assert_called_with(1)
    finally:
        os.chdir(orig_cwd)


@patch('gem2caom2.composable.GemClientCollection')
@patch('gem2caom2.pull_augmentation.retrieve_headers')
@patch('gem2caom2.gemini_metadata.retrieve_json')
@patch('gem2caom2.data_source.query_endpoint_session')
@patch('caom2pipe.client_composable.query_tap_client')
def test_run_incremental_diskfiles(
    tap_mock,
    query_mock,
    json_mock,
    header_mock,
    clients_mock,
    test_config,
    tmp_path,
    change_test_dir,
):
    query_mock.side_effect = gem_mocks.mock_query_endpoint_4
    tap_mock.side_effect = gem_mocks.mock_query_tap
    json_mock.side_effect = gem_mocks.mock_retrieve_json
    header_mock.side_effect = gem_mocks._mock_retrieve_headers_37
    clients_mock.return_value.metadata_client.read.side_effect = gem_mocks.read_mock_37

    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'testproxy.pem'
    test_config.task_types = [TaskType.INGEST]
    test_config.interval = 28800
    gem_mocks.use_test_data_cache(test_config)
    Config.write_to_file(test_config)
    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    _write_state(
        prior_timestamp='2024-10-20 20:03:00.000000',
        end_timestamp=datetime(year=2024, month=11, day=1, hour=23, minute=3, second=0),
        fqn=test_config.state_fqn,
    )

    test_result = composable._run_incremental_diskfiles()
    assert test_result is not None, 'expect result'
    assert test_result == -1, 'expect failure, no metadata mocking results set up'
    assert not clients_mock.return_value.data_client.put.called, 'data put called'
    # 16 is the number of records processed
    assert clients_mock.return_value.metadata_client.read.called, 'meta read called'
    assert clients_mock.return_value.metadata_client.read.call_count == 16, 'meta read count'
    assert not clients_mock.return_value.metadata_client.update.called, 'meta update called'
    assert not tap_mock.called, 'tap called'
    assert query_mock.called, 'query endpoint session called'
    # https://archive.gemini.edu/diskfiles/entrytimedaterange=2024-10-20T20:03:00--2024-11-01T23:03:00
    assert query_mock.call_count == 1, 'query endpoint session count'
    assert json_mock.called, 'json mock called'
    assert json_mock.call_count == 16, 'json mock count'
    assert header_mock.called, 'header mock called'
    assert header_mock.call_count == 16, 'header mock count'


@patch('gem2caom2.gemini_metadata.GeminiOrganizeExecutesRunnerMeta.do_one')
@patch('gem2caom2.composable.GemClientCollection')
@patch('gem2caom2.data_source.query_endpoint_session')
@patch('caom2pipe.client_composable.query_tap_client')
def test_run_incremental_diskfiles_limit(
    tap_mock,
    query_mock,
    clients_mock,
    do_one_mock,
    test_config,
    tmp_path,
    change_test_dir,
):

    # test that, when the limit of files is hit, the datetime for the increment is that from the time of the record
    # at the limit, not the upper end of the timebox

    query_mock.side_effect = gem_mocks.mock_query_endpoint_5
    tap_mock.side_effect = gem_mocks.mock_query_tap
    clients_mock.return_value.metadata_client.read.side_effect = gem_mocks.read_mock_37
    do_one_mock.return_value = (0, None)

    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'testproxy.pem'
    test_config.task_types = [TaskType.INGEST]
    test_config.interval = 60
    test_config.logging_level = 'WARNING'
    gem_mocks.use_test_data_cache(test_config)
    Config.write_to_file(test_config)
    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    _write_state(
        prior_timestamp='2024-08-28 17:05:00.000000',
        end_timestamp=datetime(year=2024, month=8, day=28, hour=18, minute=5, second=0),
        fqn=test_config.state_fqn,
    )

    test_result = composable._run_incremental_diskfiles()
    assert test_result is not None, 'expect result'
    assert test_result == 0, 'expect success'
    assert query_mock.called, 'query endpoint session called'
    # https://archive.gemini.edu/diskfiles/entrytimedaterange=2024-10-20T20:03:00--2024-11-01T23:03:00
    assert query_mock.call_count == 2, 'query endpoint session count'
    query_mock.assert_has_calls(
        [
            call(
                'https://archive.gemini.edu/diskfiles/NotFail/notengineering/not_site_monitoring/canonical/'
                'entrytimedaterange=2024-08-28T17:05:00--2024-08-28T18:05:00',
                ANY,
            ),
            call(
                'https://archive.gemini.edu/diskfiles/NotFail/notengineering/not_site_monitoring/canonical/'
                'entrytimedaterange=2024-08-28T17:07:32--2024-08-28T18:05:00',
                ANY,
            ),
        ]
    )
    assert not clients_mock.return_value.data_client.put.called, 'data put called'
    assert not clients_mock.return_value.metadata_client.read.called, 'meta read called'
    assert not clients_mock.return_value.metadata_client.update.called, 'meta update called'
    assert not tap_mock.called, 'tap called'

    test_state_post = State(test_config.state_fqn, zone=timezone.utc)
    assert test_state_post.get_bookmark(test_config.bookmark) == datetime(2024, 8, 28, 18, 5, 0, 0), 'saved state'


def _write_state(prior_timestamp=None, end_timestamp=None, fqn=STATE_FILE):
    # to ensure at least one spin through the execution loop, test case
    # must have a starting time greater than one config.interval prior
    # to 'now', default interval is 10 minutes
    if prior_timestamp is None:
        prior_s = datetime.now(tz=timezone.utc).timestamp() - 15 * 60
    else:
        if type(prior_timestamp) is float:
            prior_s = prior_timestamp
        else:
            prior_s = make_datetime(prior_timestamp)
    test_start_time = make_datetime(prior_s)
    if end_timestamp is None:
        test_bookmark = {
            'bookmarks': {
                'gemini_timestamp': {
                    'last_record': test_start_time,
                },
            },
        }
    else:
        assert isinstance(end_timestamp, datetime), 'end_timestamp wrong type'
        test_bookmark = {
            'bookmarks': {
                'gemini_timestamp': {
                    'last_record': test_start_time,
                    'end_timestamp': end_timestamp,
                },
            },
        }
    write_as_yaml(test_bookmark, fqn)


def _write_cert():
    fqn = f'{gem_mocks.TEST_DATA_DIR}/cadcproxy.pem'
    if not os.path.exists(fqn):
        with open(fqn, 'w') as f:
            f.write('cadc proxy content')
