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
import re

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

from gem2caom2.util import Inst


__all__ = ['FilterMetadataCache']


class FilterMetadata(object):
    # DB - 27-02-19
    # I’ve approached Gemini spectra like I have for DAO but the header
    # content for DAO spectra permits an accurate estimate of the upper/lower
    # wavelengths (a few % accuracy, at least for data from the last decade
    # or so).  It’s more a ‘ballpark’ estimate for Gemini even for simple
    # long-slit instruments like PHOENIX and NIRI.  And this approach isn’t
    # very consistent for spectrographs like GRACES or bHROS or GPI.  e.g.
    # for GRACES the WCS I’ve built  suggests that the full spectral range
    # is covered from one end of the detector to the other but in reality
    # there are perhaps 40 individual spectral orders stacked one above
    # the other each covering a few % of the full wavelength range with a
    # bit of overlap between successive orders.  I’ve never actually tried
    # to do a spectral cutout but the returned cutouts for unprocessed
    # GRACES data wouldn’t be valid.
    #
    # Just talked with Chris briefly as well.  For consistency and to avoid
    # misleading users it’s likely best to use the ‘range’ approach.
    #
    # microns for the units, not nm.  And second line would be something like:
    #
    # axis.range = CoordRange1D(RefCoord(0.5,crval1),
    #                           RefCoord(1.5,(crval1+bandpass)))
    #
    # and crval1 is the lower wavelength.

    def __init__(self, instrument=None, filter_name=None):
        self.central_wl = None
        self.bandpass = None
        self.resolving_power = None
        self.instrument = instrument
        self._filter_name = filter_name

    def __str__(self):
        return (
            f'central_wl: {self._central_wl}\n'
            f'bandpass: {self._bandpass}\n'
            f'resolving_power: {self._resolving_power}'
            f'filter_name: {self._filter_name}'
        )

    @property
    def central_wl(self):
        """Central wavelength for a filter."""
        return self._central_wl

    @central_wl.setter
    def central_wl(self, value):
        self._central_wl = value

    @property
    def bandpass(self):
        """Width of a filter."""
        return self._bandpass

    @bandpass.setter
    def bandpass(self, value):
        self._bandpass = value

    @property
    def filter_name(self):
        return self._filter_name

    @property
    def resolving_power(self):
        if self.instrument in ['Phoenix', 'TEXES']:
            return None
        elif self._resolving_power is None:
            if self.instrument == 'NIRI':
                return None
            else:
                self.adjust_resolving_power()
                return self._resolving_power
        else:
            return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        self._resolving_power = value

    def adjust_bandpass(self, variance):
        self.bandpass = (self.central_wl + variance) - (self.central_wl - variance)

    def set_bandpass(self, w_max, w_min):
        self.bandpass = w_max - w_min

    def set_central_wl(self, w_max, w_min):
        self.central_wl = (w_max + w_min) / 2.0

    def set_resolving_power(self, w_max, w_min):
        self.resolving_power = (w_max + w_min) / (2 * self.bandpass)

    def adjust_resolving_power(self):
        # the formula for direct imaging data
        self.resolving_power = self.central_wl / self.bandpass


class FilterMetadataCache:
    """Filter information that has a life-space of the application."""

    _filter_replacements = None

    def __init__(self, session):
        self._svo_session = session
        self._fm = {}
        self._logger = logging.getLogger(self.__class__.__name__)

    def __str__(self):
        msg = ''
        for key, value in self._fm.items():
            msg = f'{msg},\n{key}:{value}'
        msg = f'{msg}\n'
        return msg

    @property
    def svo_session(self):
        return self._svo_session

    @property
    def fm(self):
        return self._fm

    def filter_metadata(self, instrument, filters):
        """
        For the given instrument and filters, go to the SVO Filter Profile
        Service

        http://svo2.cab.inta-csic.es/svo/theory/fps/

        and return energy metadata.

        :param instrument: The instrument name.
        :param filters: The filter name.
        :return: FilterMetadata instance, or None, if there's no SVO information
            for the filter.
        """

        try:
            filter_names = filters.split('+')
            # use detector maximums as defaults
            wl_min = 0.0
            wl_max = 100000.0
            width_min = 100000.0
            wl_width = wl_max - wl_min
            wl_eff = (wl_max + wl_min) / 2.0

            # does the filter exist at SVO?
            filter_name_found = False

            for index in filter_names:
                filter_name = index.strip()
                filter_id = f"{instrument}.{filter_name}"
                # VERB=0 parameter means the smallest amount returned
                if instrument == 'Flamingos':
                    url = f"{ac.SVO_URL}KPNO/{filter_id}&VERB=0"
                else:
                    url = f"{ac.SVO_URL}Gemini/{filter_id}&VERB=0"

                # Open the URL and fetch the VOTable document.
                # Some Gemini filters in SVO filter database have bandpass info
                # only for 'w'arm filters.  First check for filter without 'w'
                # appended to the ID (which I assume means bandpass is for cold
                # filter), then search for 'w' if nothing is found...
                votable, error_message = ac.get_vo_table_session(url, self._svo_session)
                if not votable:
                    if instrument == 'Flamingos':
                        url = f"{ac.SVO_URL}KPNO/{filter_id}w&VERB=0"
                    else:
                        url = f"{ac.SVO_URL}Gemini/{filter_id}w&VERB=0"
                    votable, error_message = ac.get_vo_table_session(url, self._svo_session)
                if not votable:
                    logging.error(f'Unable to download SVO filter data from {url} because {error_message}')
                    continue

                # DB - 14-04-19 After discussion with a few others use the
                # wavelength lookup values “WavelengthCen” and “FWHM” returned
                # from the SVO. Looking at some of the IR filters
                # use the more common “WavelengthCen” and “FWHM” values that
                # the service offers.

                filter_name_found = True
                wl_width = votable.get_field_by_id('FWHM').value
                wl_eff = votable.get_field_by_id('WavelengthCen').value
                w_min = wl_eff - wl_width / 2.0
                w_max = wl_eff + wl_width / 2.0

                if w_min > wl_min:
                    wl_min = w_min
                if w_max < wl_max:
                    wl_max = w_max
                if wl_width < width_min:
                    width_min = wl_width

            if filter_name_found:
                local_fm = FilterMetadata(instrument, filter_names)
                # SVO filter units are angstroms, Gemini CAOM2 spectral wcs is
                # microns
                local_fm.central_wl = wl_eff / 1.0e4
                local_fm.bandpass = wl_width / 1.0e4
                logging.info(f'Filter(s): {filter_names}  MD: {local_fm.central_wl}, {local_fm.bandpass}')
                return local_fm
            else:
                return None
        except Exception as e:
            logging.error(e)
            import traceback
            tb = traceback.format_exc()
            logging.debug(tb)

    def get_filter_metadata(self, instrument, filter_name, telescope):
        """A way to lazily initialize all the filter metadata reads from SVO."""
        self._logger.debug(
            f'Begin get_filter_metadata with instrument {instrument} '
            f'filter name {filter_name} telescope {telescope}'
        )
        repaired_inst = FilterMetadataCache._repair_instrument_name_for_svo(instrument, telescope)
        repaired_filters = FilterMetadataCache._repair_filter_name_for_svo(instrument, filter_name)
        self._logger.debug(f'Find information for filter {repaired_filters} on instrument {repaired_inst}')
        if repaired_filters is None:
            # nothing to look up, try something else
            return None
        if repaired_inst in self._fm and repaired_filters in self._fm[repaired_inst]:
            result = self._fm[repaired_inst][repaired_filters]
            if result is not None:
                result.adjust_resolving_power()
        else:
            result = self.filter_metadata(repaired_inst, repaired_filters)
            if result is None:
                self._logger.warning(f'No filter information found for {repaired_inst}:{repaired_filters}')

            # keep the None result even if there is no filter information to stop repeated lookups
            if repaired_inst in self._fm:
                temp = self._fm[repaired_inst]
                temp[repaired_filters] = result
            else:
                self._fm[repaired_inst] = {repaired_filters: result}
        return result

    @staticmethod
    def _filter_replacement_map():
        """Load filter_replacement from the pipeline cache on first use."""
        if FilterMetadataCache._filter_replacements is None:
            replacements = mc.Cache().get_from('filter_replacement')
            if not isinstance(replacements, dict):
                raise mc.CadcException('filter_replacement must be a map of instrument name to filter name.')
            legal = {member.value for member in Inst if member.value is not None}
            for key in replacements:
                if key not in legal:
                    raise mc.CadcException(f'Unknown filter_replacement instrument {key}.')
            FilterMetadataCache._filter_replacements = replacements
        return FilterMetadataCache._filter_replacements

    @staticmethod
    def _repair_filter_name_for_svo(instrument, filter_names):
        """
        Filter names from JSON/headers are not necessarily the same
        as the filter names used by the SVO Filter service. Correlate
        the two here.

        DB - 02-04-19 - strip the bar code from the filter names

        :param instrument what repairs to apply
        :param filter_names the Gemini version, which may include multiple names
            separated by '+'
        :return filter_name the SVO version
        """
        replacements = FilterMetadataCache._filter_replacement_map()
        table = None
        if isinstance(instrument, Inst) and instrument.value is not None:
            table = replacements.get(instrument.value)

        result = []
        for filter_name in filter_names.split('+'):
            temp = filter_name
            if instrument is Inst.NIRI:
                temp = re.sub(r'con', 'cont', temp)
                temp = re.sub(r'_', '-', temp)
                temp = re.sub('\\(', '', temp)
                temp = re.sub('\\)', '', temp)
                temp = temp if table is None else table.get(temp, temp)
            elif instrument in [Inst.TRECS, Inst.MICHELLE]:
                temp = filter_name.split('-')
                if len(temp) > 0:
                    temp = temp[0]
                temp = temp if table is None else table.get(temp, temp)
            elif instrument is Inst.HRWFS:
                # “ND” in the filter name means ‘neutral density’.  Ignore any
                # of these as they have no impact on the transmitted wavelengths
                # - I think #159 was the only one delivered according to
                # http://www.gemini.edu/sciops/telescope/acqcam/acqFilterList.html.
                # Acqcam/hrwfs was used mainly to look for rapid variability in
                # bright, stellar objects that were really too bright for an 8'
                # telescope and would have saturated the detector without an ND
                # filter.
                if temp.startswith('ND'):
                    continue
                temp = temp[0]
            elif instrument is Inst.F2:
                if temp.startswith('J-lo'):
                    temp = 'Jlow'
                temp = temp if table is None else table.get(temp, temp)
            else:
                temp = temp if table is None else table.get(temp, temp)
            if temp is not None:
                result.append(temp)
        if len(result) > 0:
            return '+'.join(i for i in result)
        else:
            return None

    @staticmethod
    def _repair_instrument_name_for_svo(instrument, telescope=None):
        """
        Instrument names from JSON/headers are not necessarily the same
        as the instrument names used by the SVO Filter service. Correlate
        the two here.
        :param instrument the Gemini version
        :return instrument the SVO version
        """
        result = instrument.value
        if instrument is Inst.HRWFS:
            if telescope is None:
                raise mc.CadcException(f'{instrument}: No observatory information.')
            else:
                if 'Gemini-South' == telescope:
                    result = 'AcqCam-S'
                else:
                    result = 'AcqCam-N'
        elif instrument is Inst.F2:
            result = 'Flamingos2'
        elif instrument is Inst.FLAMINGOS:
            result = 'Flamingos'
        return result
