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

import os
import warnings

import pytest

from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import FITSFixedWarning
from caom2utils import data_util
from caom2pipe.manage_composable import CadcException, ExecutionReporter2, read_obs_from_file, TaskType
from gem2caom2.util import Inst
from gem2caom2 import fits2caom2_augmentation, gemini_metadata, obs_file_relationship, pull_augmentation, svofps
from gem2caom2.gem_name import GemName
from gem2caom2.program_metadata import MDContext, PIMetadata

from unittest.mock import ANY, Mock, patch
import gem_mocks


pytest.main(args=['-s', os.path.abspath(__file__)])


def pytest_generate_tests(metafunc):
    if os.path.exists(gem_mocks.TEST_DATA_DIR):

        file_list = []
        for ii in [
            Inst.GMOS,
            Inst.NIRI,
            Inst.GPI,
            Inst.F2,
            Inst.GSAOI,
            Inst.NICI,
            Inst.TRECS,
            Inst.MICHELLE,
            Inst.GRACES,
            Inst.NIFS,
            Inst.GNIRS,
            Inst.PHOENIX,
            Inst.FLAMINGOS,
            Inst.HRWFS,
            Inst.HOKUPAA,
            Inst.OSCIR,
            Inst.BHROS,
            Inst.CIRPASS,
            Inst.TEXES,
            'processed',
            Inst.ALOPEKE,
            Inst.ZORRO,
            Inst.IGRINS,
            Inst.IGRINS2,
            Inst.GHOST,
            Inst.MAROONX,
        ]:
            walk_dir = _get_inst_name(ii)
            for root, dirs, files in os.walk(f'{gem_mocks.TEST_DATA_DIR}/{walk_dir}'):
                for file in files:
                    if file.endswith(".header"):
                        file_list.append(os.path.join(root, file))

        metafunc.parametrize('test_name', file_list)


@patch('gem2caom2.pull_augmentation.PullVisitor.look_pull_and_put')
@patch('gem2caom2.program_metadata.mc.query_endpoint_session')
@patch('caom2utils.data_util.get_file_type')
@patch('gem2caom2.gemini_metadata.retrieve_headers')
@patch('gem2caom2.gemini_metadata.retrieve_json')
@patch('caom2pipe.astro_composable.get_vo_table_session')
@patch('gem2caom2.main_app.ProvenanceFinder')
def test_visitor(
    pf_mock,
    svofps_mock,
    json_mock,
    header_mock,
    file_type_mock,
    pi_mock,
    pull_mock,
    test_name,
    test_config,
    tmp_path,
    change_test_dir,
):
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    warnings.simplefilter('ignore', AstropyWarning)
    test_file_id = obs_file_relationship.remove_extensions(os.path.basename(test_name))
    svofps_mock.side_effect = gem_mocks.mock_get_votable
    pf_mock.return_value.get.side_effect = gem_mocks.mock_get_data_label
    json_mock.side_effect = gem_mocks.mock_get_obs_metadata
    file_type_mock.return_value = 'application/fits'
    pi_mock.side_effect = gem_mocks.mock_get_pi_metadata

    test_config.task_types = [TaskType.SCRAPE]
    test_config.use_local_files = True
    test_config.log_to_file = True
    test_config.data_sources = [tmp_path.as_posix()]
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.logging_level = 'ERROR'
    test_config.write_to_file(test_config)
    gem_mocks.set_logging(test_config)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    expected_fqn = f'{os.path.dirname(test_name)}/{test_file_id}.expected.xml'
    in_fqn = expected_fqn.replace('.expected.', '.in.')

    clients_mock = Mock()

    def _mock_repo_read(collection, obs_id):
        input_observation = None
        if os.path.exists(in_fqn):
            input_observation = read_obs_from_file(in_fqn)
        return input_observation

    clients_mock.metadata_client.read.side_effect = _mock_repo_read

    def _read_header_mock(ignore1, ignore2, ignore3, ignore4, ignore5):
        return data_util.get_local_file_headers(test_name)

    header_mock.side_effect = _read_header_mock

    actual_fqn = expected_fqn.replace('expected', 'actual')
    if os.path.exists(actual_fqn):
        os.unlink(actual_fqn)

    test_reporter = ExecutionReporter2(test_config)
    filter_cache = svofps.FilterMetadataCache(svofps_mock)
    pi_metadata = PIMetadata(gemini_session=Mock())
    md_context = MDContext(filter_cache, pi_metadata)
    test_subject = gemini_metadata.GeminiMetaVisitRunnerMeta(
        clients_mock, test_config, [pull_augmentation, fits2caom2_augmentation], test_reporter
    )
    storage_name = GemName(test_name, md_context)
    if gem_mocks.LOOKUP[test_file_id][1] == Inst.ZORRO:
        storage_name.obs_id = test_file_id[:-1]
    else:
        storage_name.obs_id = gem_mocks.LOOKUP[test_file_id][0]
    context = {'storage_name': storage_name}
    try:
        test_subject.execute(context)
    except CadcException as e:
        if storage_name.file_name in ['N20220915S0113.fits', 'S20230301S0170.fits']:
            assert test_reporter._observable.rejected.is_mystery_value(
                storage_name.file_name
            ), 'expect rejected mystery value record'
        else:
            raise e

    gem_mocks.compare(expected_fqn, actual_fqn, test_subject._observation)

    if test_subject._observation.observation_id in ['GS-2022B-Q-235-137-045', 'GS-2023A-LP-103-23-017']:
        assert test_reporter._observable.rejected.is_bad_metadata(storage_name.file_name), 'expect rejected record'
    else:
        assert not test_reporter._observable.rejected.is_bad_metadata(
            storage_name.file_name
        ), 'expect no rejected record'

    assert json_mock.called, 'json mock'
    json_mock.assert_called_with(storage_name.file_id, ANY, ANY), 'json mock args'
    if storage_name.file_name in [
        'N20120825S0597_arc.fits',
        '2007sep15_0001.fits',
        'TX20131117_flt.3002.fits',
        'TX20131117_raw.3002.fits',
    ]:
        # the proprietary files
        pull_mock.assert_not_called(), 'pull mock not called'
    else:
        pull_mock.assert_called_once(), 'pull mock called once'

    # assert False


def _get_inst_name(inst):
    walk_dir = inst
    if inst != 'processed' and isinstance(inst, Inst):
        walk_dir = inst.value
    return walk_dir
