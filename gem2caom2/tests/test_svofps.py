# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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
#  : 4 $
#
# ***********************************************************************
#

from gem2caom2.util import Inst
from gem2caom2 import svofps

from caom2pipe.manage_composable import CadcException
from mock import patch, Mock, ANY
import gem_mocks
import pytest


@pytest.fixture(autouse=True)
def _chdir_test_data(monkeypatch, test_data_dir):
    monkeypatch.chdir(test_data_dir)


test_subjects = {
    'H2v=2-1S1_G0220': [Inst.NIRI, 'H2S1v2-1-G0220'],
    'Kprime_G0206': [Inst.NIRI, 'Kprime-G0206'],
    'H2Oice204_G0242': [Inst.NIRI, 'H2Oice2045-G0242'],
    'Jcon(121)_G0232': [Inst.NIRI, 'Jcont1207-G0232'],
    'Bra_G0238': [Inst.NIRI, 'BrAlpha-G0238'],
    'Bracont_G0237': [Inst.NIRI, 'BrAlphaCont-G0237'],
    'Brgamma_G0218': [Inst.NIRI, 'BrG-G0218'],
    'Jcon1065_G0239': [Inst.NIRI, 'Jcont1065-G0239'],
    'hydrocarb_G0231': [Inst.NIRI, 'hydrocarbon-G0231'],
}


def test_repair_filter_name():
    for ii in test_subjects:
        test_result = svofps.FilterMetadataCache._repair_filter_name_for_svo(test_subjects[ii][0], ii)
        assert test_result == test_subjects[ii][1], 'wrong value'


def test_repair_filter_name_miss_policy():
    # NICI and Alopeke/Zorro drop a name that has no replacement.
    # Alopeke and Zorro share one map.
    cases = [
        (Inst.NICI, 'J', 'J'),
        (Inst.ALOPEKE, 'open', 'open'),
        (Inst.ALOPEKE, 'Red-832', 'EO_832'),
        (Inst.ZORRO, 'Red-832', 'EO_832'),
    ]
    for instrument, filter_name, expected in cases:
        test_result = svofps.FilterMetadataCache._repair_filter_name_for_svo(instrument, filter_name)
        assert test_result == expected, f'{instrument} {filter_name}'


@patch('caom2pipe.astro_composable.get_vo_table_session')
def test_get_filter_metadata(get_vo_mock):
    get_vo_mock.side_effect = gem_mocks.mock_get_votable
    cache = Mock()
    test_subject = svofps.FilterMetadataCache(Mock())
    test_result = test_subject.get_filter_metadata(Inst.NIRI, 'filters', 'telescope_name')
    assert get_vo_mock.call_count == 2, 'wrong number of calls'
    get_vo_mock.assert_called_with(
        'http://svo2.cab.inta-csic.es/svo/theory/fps3/fps.php?ID=Gemini/' 'NIRI.filtersw&VERB=0',
        ANY,
    ), 'wrong call args'
    assert test_result is None, 'do not expect a result'
    # do the same thing again, check that the result has been cached
    test_result = test_subject.get_filter_metadata(Inst.NIRI, 'filters', 'telescope_name')
    assert get_vo_mock.call_count == 2, 'wrong number of calls'
    assert test_result is None, 'do not expect a result this time either'


def test_unknown_instrument_key(monkeypatch, tmp_path):
    svofps.FilterMetadataCache._filter_replacements = None
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'cache.yml').write_text('filter_replacement:\n  NotAnInstrument:\n    A: B\n')
    (tmp_path / 'config.yml').write_text(f"working_directory: '{tmp_path}'\ncache_file_name: cache.yml\n")
    with pytest.raises(CadcException, match='NotAnInstrument'):
        svofps.FilterMetadataCache._repair_filter_name_for_svo(Inst.NIRI, 'Kprime_G0206')
