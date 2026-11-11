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
"""
Notes on the GEM archive/GEMINI collection:

# TODO - none of this is really true anymore, so figure out what needs
to be captured here for documentation ..... ;)

1. Must use the file name as the starting point for work, because that's
what is coming back from the ad query, and the ad query is what is being
used to trigger the work.

2. Must find the observation ID value from the file header information,
because the observation ID is how to get the existing CAOM instance.

3. Artifact URIs in existing observations reference the gemini schema, not
the ad schema.

4. TODO - what happens to obtaining the observation ID, if a preview
file is retrieved prior to the FITS file being retrieved?

Because of this, make the GemName a class that, standing on it's own,
can retrieve the observation ID value from the headers for a file.

"""

import logging
import math
import re
import traceback

from astropy.time import Time

from caom2 import CalibrationLevel, Chunk, ProductType
from caom2 import TypedList, DerivedObservation, DataProductType
from caom2 import ObservationIntentType, TargetType, CoordAxis1D, Axis
from caom2 import SpectralWCS, RefCoord, CoordRange1D
from caom2utils.caom2blueprint import update_artifact_meta
from caom2utils.blueprints import _to_float
from caom2utils.wcs_parsers import FitsWcsParser
from caom2utils.data_util import get_file_type
from caom2pipe import manage_composable as mc
from caom2pipe import caom_composable as cc
from caom2pipe import astro_composable as ac

import gem2caom2.obs_file_relationship as ofr
from gem2caom2.gem_name import GemName
from gem2caom2 import program_metadata, svofps
from gem2caom2.gemini_metadata import GeminiMetadataLookupStorageName, ProvenanceFinder
from gem2caom2.util import Inst


__all__ = ['GeminiMapping', 'mapping_factory']

# DB 01-04-21
# PINHOLE is CALIBRATION
# DB 03-06-21
# OBSCLASS = dayCal datasets should all have an intent of calibration.
CAL_VALUES = [
    'GCALflat',
    'Bias',
    'BIAS',
    'Twilight',
    'Ar',
    'FLAT',
    'flat',
    'ARC',
    'Domeflat',
    'DARK',
    'dark',
    'gcal',
    'ZERO',
    'SLIT',
    'slit',
    'PINHOLE',
    'dayCal',
]


# DB - 04-02-19 - strip out anything with 'pupil' as it doesn't affect
# energy transmission
#
# DB 24-04-19
# Other instruments occasionally have ND filters in the beam.
FILTER_VALUES_TO_IGNORE = ['open', 'invalid', 'pupil', 'clear', 'nd']

# GPI radius == 2.8 arcseconds, according to Christian Marois via DB 02-07-19
# DB - 18-02-19 - Replace “5.0” for “2.8" for GPI field of view.
# NOTE:  To be more accurate for GRACES this size could be reduced
# from 5" to 1.2" since that’s the size of the fibre.
# OSCIR - http://www.gemini.edu/sciops/instruments/oscir/oscirIndex.html
# bHROS - DB - 20-02-19 - bHROS ‘bounding box’ is only 0.9".
#                         A very small fibre.
# HOKUPAA - http://www.gemini.edu/sciops/instruments/uhaos/uhaosIndex.html
# NIFS - DB - 04-03-19 - hard-code 3" FOV
# TEXES - DB - 07-03-19 - cd11=cd22 = 5.0/3600.0
# IGRINS - DB - 21-07-21 - hard-code 5" FOV
# GHOST - DB - 12-01-24 - hard-code to sqrt(0.94)"
# MAROONX - https://gemini.edu/instrumentation/maroon-x/capability
RADIUS_LOOKUP = {
    Inst.GPI: 2.8 / 3600.0,  # units are arcseconds
    Inst.GRACES: 1.2 / 3600.0,
    Inst.PHOENIX: 5.0 / 3600.0,
    Inst.OSCIR: 0.0890 / 3600.0,
    Inst.HOKUPAA: 4.0 / 3600.0,
    Inst.BHROS: 0.9 / 3600.0,
    Inst.NIFS: 3.0 / 3600.0,
    Inst.TEXES: 5.0 / 3600.0,
    Inst.IGRINS: 5.0 / 3600.0,
    Inst.GHOST: math.sqrt(0.94) / 3600.0,
    Inst.MAROONX: math.sqrt(0.77) / 3600.0,
}


class GeminiValueRepair(mc.ValueRepairCache):

    VALUE_REPAIR = {
        'observation.type': {
            'dark': 'DARK',
            'Dark': 'DARK',
            'object': 'OBJECT',
            'flat': 'FLAT',
            'Flat': 'FLAT',
        },
        'chunk.position.axis.axis1.ctype': {
            'RA--TAN': 'RA---TAN',
        },
        # WF - 02-02-23
        # APPT is a typo, repair to GAPPT
        'chunk.position.coordsys': {'APPT': 'GAPPT'},
    }

    def __init__(self):
        self._value_repair = GeminiValueRepair.VALUE_REPAIR
        self._key = None
        self._values = None
        self._logger = logging.getLogger(self.__class__.__name__)


class GeminiMapping(cc.TelescopeMapping2):

    # value repair cache
    value_repair = GeminiValueRepair()

    def __init__(
        self,
        storage_name,
        lookup,
        instrument,
        clients,
        reporter,
        observation,
        config,
        provenance_finder,
    ):
        super().__init__(storage_name, clients, reporter, observation, config)
        self._instrument = instrument
        self._lookup = lookup
        self._filter_cache = self._storage_name.md_context.filter_cache
        self._provenance_finder = provenance_finder
        self.fm = None

    def get_art_content_checksum(self, ext):
        return f'md5:{self._lookup.data_md5(self._storage_name.file_uri)}'

    def get_art_content_length(self, ext):
        return self._lookup.data_size(self._storage_name.file_uri)

    def get_art_content_type(self, ext):
        return get_file_type(self._storage_name.file_name)

    def get_calibration_level(self, ext):
        result = CalibrationLevel.RAW_STANDARD
        reduction = self._lookup.reduction(self._storage_name.file_uri)
        if reduction is not None and ('PROCESSED' in reduction or 'PREPARED' in reduction):
            result = CalibrationLevel.CALIBRATED
        return result

    def get_cd11(self, ext, keyword='CDELT1'):
        result = self._headers[ext].get(keyword)
        if result is None:
            result = RADIUS_LOOKUP.get(self._instrument)
        return result

    def get_cd22(self, ext):
        return self.get_cd11(ext, 'CDELT2')

    def get_crpix1(self, ext):
        return 1.0

    def get_crpix2(self, ext):
        return 1.0

    def get_art_product_type(self, ext):
        """
        If obsclass is unknown then CAOM2 ProductType is set to CALIBRATION.
        This should only effect early data when OBSCLASS was not in the JSON
        summary metadata or FITS headers.
        """
        data_label = self._lookup.data_label(self._storage_name.file_uri)
        obs_type = self._lookup.observation_type(self._storage_name.file_uri)
        obs_class = self._lookup.observation_class(self._storage_name.file_uri)

        self._logger.debug(f'obs type is {obs_type} obs class is {obs_class} for {data_label}')
        if obs_type is not None and obs_type == 'MASK':
            result = ProductType.AUXILIARY
        elif obs_class is None:
            # 'unknown' is the value of observation_class for CIRPASS
            if obs_type is not None and obs_type in ['OBJECT', 'unknown']:
                if data_label is not None and 'CAL' in data_label:
                    result = ProductType.CALIBRATION
                else:
                    result = ProductType.SCIENCE
            else:
                result = ProductType.CALIBRATION
        elif obs_class is not None and obs_class == 'science':
            result = ProductType.SCIENCE
        else:
            result = ProductType.CALIBRATION
        return result

    def get_data_product_type(self, ext):
        mode = self._lookup.mode(self._storage_name.file_uri)
        obs_type = self._lookup.observation_type(self._storage_name.file_uri)
        if mode is None:
            raise mc.CadcException(
                f'{self._instrument}: No mode information found for {self._storage_name.file_name}'
            )
        elif (mode == 'imaging') or (obs_type is not None and obs_type == 'MASK'):
            result = DataProductType.IMAGE
        else:
            result = DataProductType.SPECTRUM
        return result

    def get_data_release(self, ext):
        """
        Determine the plane-level data release date.
        :return: The Plane release date, or None if not found.
        """
        # every instrument has a 'release' keyword in the JSON summary not every instrument (Michelle) has a
        # RELEASE keyword in the appropriate headers
        return self._lookup.release(self._storage_name.file_uri)

    def get_dec(self, ext):
        return self._lookup.dec(self._storage_name.file_uri)

    def get_exposure(self, ext):
        """
        Calculate the exposure time. EXPTIME in the header is not always the
        total exposure time for some of the IR instruments.  For these EXPTIME
        is usually the individual exposure time for each chop/nod and a bunch
        of chops/nods are executed for one observation.  Gemini correctly
        allows for this in the json ‘exposure_time’ value.
        """
        return self._lookup.exposure_time(self._storage_name.file_uri)

    def get_meta_release(self, ext):
        """
        Determine the metadata release date (Observation and Plane-level).
        :return: The Observation/Plane release date, or None if not found.
        """
        # TODO - this is a bit different ROTFL
        if self._headers:
            # DB 21-08-19
            # If PROP_MD is T, use JSON ‘release’ value for metadata release
            # date. If no PROP_MD present or value is F use the JSON
            # ut_datetime value.
            prop_md = self._headers[ext].get('PROP_MD')
            meta_release = None
            if prop_md is None or prop_md is False or prop_md == 'F':
                meta_release = self._lookup.ut_datetime(self._storage_name.file_uri)
            if meta_release is None or prop_md and prop_md == 'T' or prop_md:
                # a late value is better than None, since that makes the Observation unfindable
                meta_release = self._lookup.release(self._storage_name.file_uri)
        else:
            # GenericParser, so no headers retrieved from archive.gemini.edu,
            # probably a 403 being returned by the site, assume proprietary
            meta_release = self._lookup.release(self._storage_name.file_uri)

        return mc.make_datetime(meta_release)

    def get_obs_intent(self, ext):
        result = ObservationIntentType.CALIBRATION
        data_label = self._lookup.data_label(self._storage_name.file_uri)
        obs_class = self._lookup.observation_class(self._storage_name.file_uri)
        self._logger.debug(f'observation_class is {obs_class} for {data_label}')
        if obs_class is None:
            obs_type = self._lookup.observation_type(self._storage_name.file_uri)

            self._logger.debug(f'observation_type is {obs_type} for {data_label}')
            if obs_type is None:
                self._logger.debug(f'data_label is {data_label}')
                if (
                    data_label is None
                    or (data_label is not None and '-CAL' not in data_label)
                    and self._headers[ext] is not None
                ):
                    object_value = self._headers[ext].get('OBJECT')
                    self._logger.debug(f'object_value is {object_value} for {data_label}')
                    if object_value is not None:
                        if object_value in CAL_VALUES:
                            result = ObservationIntentType.CALIBRATION
                        else:
                            # check that an individual cal value is not part
                            # of the object_value
                            cal_value_found = False
                            for ii in CAL_VALUES:
                                if ii in object_value:
                                    result = ObservationIntentType.CALIBRATION
                                    cal_value_found = True
                                    break
                            if not cal_value_found:
                                result = ObservationIntentType.SCIENCE
                else:
                    if '-CAL' in data_label:
                        result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
            else:
                if obs_type in CAL_VALUES:
                    result = ObservationIntentType.CALIBRATION
                else:
                    result = ObservationIntentType.SCIENCE
        elif 'science' in obs_class:
            result = ObservationIntentType.SCIENCE
        return result

    def get_obs_type(self, ext):
        result = self._lookup.observation_type(self._storage_name.file_uri)
        obs_class = self._lookup.observation_class(self._storage_name.file_uri)
        if obs_class is not None and (obs_class == 'acq' or obs_class == 'acqCal'):
            result = 'ACQUISITION'
        return result

    def get_proposal_id(self, ext):
        return self._lookup.program_id(self._storage_name.file_uri)

    def get_provenance_keywords(self, ext):
        """
        DB https://github.com/opencadc-metadata-curation/gem2caom2/issues/12
        Currently there is no CAOM2 metadata that enables a user to distinguish
        different spectroscopic modes of GMOS-N/S data. e.g. long-slit vs. IFU
        vs. MOS (multi-object spectroscopy).

        To enable this with at least a TAP query it would be useful to modify the
        gem2caom2 code to add an Instrument.keywords value for GMOS-N/S spectra.

        The jsonsummary 'mode' value is likely sufficient for this. It is
        supposed to provide values of "imaging, spectroscopy, LS (Longslit
        Spectroscopy), MOS (Multi Object Spectroscopy) or IFS (Integral Field
        Spectroscopy)". There is likely no reason to change these values but
        simply use them for the value of Instrument.keywords.
        """
        return self._lookup.mode(self._storage_name.file_uri)

    def get_provenance_last_executed(self, ext):
        def breakout(comments):
            result = None
            temp = comments.split('\n')
            if len(temp) > 6 and 'HST' in temp[6]:
                result = mc.make_datetime(temp[6])
            return result

        return self._get_provenance_breakout(ext, breakout)

    def get_provenance_producer(self, ext):
        cal_level = self.get_calibration_level(ext)
        if cal_level in [
            CalibrationLevel.CALIBRATED,
            CalibrationLevel.PRODUCT,
            CalibrationLevel.ANALYSIS_PRODUCT,
        ]:
            comments = str(self._headers[ext].get('COMMENT'))
            result = comments.split('Processed by the')[1].split('|')[0]
        else:
            result = self._headers[ext].get('IMAGESWV')
        return result

    def get_provenance_reference(self, ext):
        def breakout(comments):
            return 'https://www.gemini.edu/instrumentation/graces/data-reduction'

        return self._get_provenance_breakout(ext, breakout)

    def get_provenance_version(self, ext):
        def breakout(comments):
            temp = comments.split('opera-')[1].split('build date')[0]
            return f'opera-{temp}'

        return self._get_provenance_breakout(ext, breakout)

    def _get_provenance_breakout(self, ext, fn):
        result = None
        cal_level = self.get_calibration_level(ext)
        if cal_level in [
            CalibrationLevel.CALIBRATED,
            CalibrationLevel.PRODUCT,
            CalibrationLevel.ANALYSIS_PRODUCT,
        ]:
            comments = str(self._headers[ext].get('COMMENT'))
            result = fn(comments)
        return result

    def get_ra(self, header):
        """
        Get the right ascension. Rely on the JSON metadata, because it's all in
        the same units (degrees).

        :param header:  The FITS header for the current extension.
        :return: ra, or None if not found.
        """
        return self._lookup.ra(self._storage_name.file_uri)

    def get_sky_coord(self, ext, ra_key, dec_key):
        ra_hours = self._headers[ext].get(ra_key)
        dec_hours = self._headers[ext].get(dec_key)
        if (
            ra_hours is None
            or dec_hours is None
            or ra_hours == 'INDEF'
            or dec_hours == 'INDEF'
            or ra_hours == 'Unknown'
            or dec_hours == 'Unknown'
        ):
            ra_deg = None
            dec_deg = None
        else:
            ra_deg, dec_deg = ac.build_ra_dec_as_deg(ra_hours, dec_hours)
        return ra_deg, dec_deg

    def get_target_moving(self, ext):
        """
        Calculate whether the Target moving.
        Non-sidereal tracking -> setting moving target to "True"

        DB, 01-08-19
        Many calibration observations are acquired with the telescope parked
        and hence not tracking at sidereal rate.
        """
        types = self._lookup.types(self._storage_name.file_uri)
        if 'NON_SIDEREAL' in types:
            return True
        else:
            return None

    def get_target_type(self, ext):
        """
        Calculate the Target TargetType
        """
        result = TargetType.FIELD
        spectroscopy = self._lookup.spectroscopy(self._storage_name.file_uri)
        if spectroscopy:
            result = TargetType.OBJECT
        return result

    def get_time_delta(self, ext):
        result = None
        exptime = self.get_exposure(ext)
        if exptime is not None:
            result = mc.convert_to_days(mc.to_float(exptime))
        return result

    def get_time_function_val(self, ext):
        """
        Calculate the Chunk Time WCS function value, in 'mjd'.
        """
        time_string = self._lookup.ut_datetime(self._storage_name.file_uri)
        result = ac.get_datetime_mjd(mc.make_datetime(time_string))
        if result and isinstance(result, Time):
            result = result.value
        return result

    def _accumulate_chunk_position_axes_blueprint(self, bp, extension):
        bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN', extension)
        bp.set('Chunk.position.axis.axis1.cunit', 'deg', extension)
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN', extension)
        bp.set('Chunk.position.axis.axis2.cunit', 'deg', extension)
        bp.set('Chunk.position.axis.function.cd11', 'get_cd11()', extension)
        bp.set('Chunk.position.axis.function.cd12', 0.0, extension)
        bp.set('Chunk.position.axis.function.cd21', 0.0, extension)
        bp.set('Chunk.position.axis.function.cd22', 'get_cd22()', extension)
        bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
        bp.set('Chunk.position.axis.function.dimension.naxis2', 1)
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix', 'get_crpix1()', extension)
        bp.set('Chunk.position.axis.function.refCoord.coord1.val', 'get_ra()', extension)
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 'get_crpix2()', extension)
        bp.set('Chunk.position.axis.function.refCoord.coord2.val', 'get_dec()', extension)
        bp.add_attribute('Chunk.position.equinox', 'EQUINOX', extension)

    def _accumulate_chunk_time_axis_blueprint(self, bp, axis):
        bp.configure_time_axis(axis)
        # The Chunk time metadata is calculated using keywords from the primary header, and the only way I could
        # figure out to access keywords in the primary is through a function.
        bp.set('Chunk.time.resolution', 'get_exposure()')
        bp.set('Chunk.time.exposure', 'get_exposure()')
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.error.syser', '1e-07')
        bp.set('Chunk.time.axis.error.rnder', '1e-07')
        bp.set('Chunk.time.axis.function.naxis', '1')
        bp.set('Chunk.time.axis.function.delta', 'get_time_delta()')
        bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
        bp.set('Chunk.time.axis.function.refCoord.val', 'get_time_function_val()')

    def _accumulate_obs_plane_artifact_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        self._logger.debug(f'Begin accumulate_blueprint for {self._storage_name.file_id}.')
        bp.set('Observation.type', 'get_obs_type()')
        bp.set('Observation.intent', 'get_obs_intent()')
        bp.set('Observation.metaRelease', 'get_meta_release()')
        bp.set('Observation.target.type', 'get_target_type()')
        bp.set('Observation.target.moving', 'get_target_moving()')
        bp.set('Observation.proposal.id', 'get_proposal_id()')

        bp.clear('Observation.algorithm.name')
        bp.set('Observation.instrument.name', self._instrument.value)
        telescope = self._lookup.telescope(self._storage_name.file_uri)
        if telescope is not None:
            if telescope is not None and 'North' in telescope:
                x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
            else:
                x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
            bp.set('Observation.telescope.geoLocationX', x)
            bp.set('Observation.telescope.geoLocationY', y)
            bp.set('Observation.telescope.geoLocationZ', z)

        bp.set('Plane.productID', self._storage_name.file_id)
        bp.set('Plane.dataProductType', 'get_data_product_type()')
        bp.set('Plane.calibrationLevel', 'get_calibration_level()')
        bp.set('Plane.metaRelease', 'get_meta_release()')
        bp.set('Plane.dataRelease', 'get_data_release()')

        bp.set('Plane.provenance.name', 'Gemini Observatory Data')
        bp.set('Plane.provenance.project', 'Gemini Archive')
        # Add IMAGESWV for GRACES
        bp.add_attribute('Plane.provenance.producer', 'IMAGESWV')
        bp.set_default('Plane.provenance.producer', 'Gemini Observatory')
        # want the un-repaired data label value for resolving at archive.gemini.edu
        data_label = self._storage_name.json_metadata.get(self._storage_name.file_uri).get('data_label')
        bp.set('Plane.provenance.reference', f'http://archive.gemini.edu/searchform/{data_label}')

        bp.set('Artifact.productType', 'get_art_product_type()')
        bp.set('Artifact.contentChecksum', 'get_art_content_checksum()')
        bp.set('Artifact.contentLength', 'get_art_content_length()')
        bp.set('Artifact.contentType', 'get_art_content_type()')
        # always see the metadata, see the data only when it's public
        bp.set('Artifact.releaseType', 'data')
        bp.set('Artifact.uri', self._storage_name.file_uri)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug(f'Begin accumulate_blueprint for {self._storage_name.file_id}.')
        self._accumulate_obs_plane_artifact_blueprint(bp)
        self._accumulate_chunk_time_axis_blueprint(bp, 3)
        self._logger.debug('Done accumulate_blueprint.')

    def update(self):

        if self._headers is None:
            # proprietary header metadata at archive.gemini.edu, so do nothing
            return self._observation

        # processed files
        if cc.is_composite(self._headers) and not isinstance(self._observation, DerivedObservation):
            self._observation = self._update_composite(self._observation)

        if self._instrument in [Inst.MICHELLE, Inst.GNIRS]:
            # CSD-2511
            # DB 16-04-19
            # The more important issue with this and other files is that they
            # contain no image extensions.  The file is downloadable from
            # the Gemini archive but their only content is the primary
            # header.   There is no pixel data.  Test for the existence of a
            # FITS extension and skip processing of a michelle file if
            # there isn’t one

            # DB 18-04-19
            #
            # For the last GNIRS file (NAXIS=0)  skip the file if it doesn’t
            # have an extension.

            if len(self._headers) == 1:
                self._logger.warning(
                    f'{self._instrument}: no image data for {self._observation.observation_id}. Cannot build an '
                    'observation.'
                )
                return None

        try:
            for plane in self._observation.planes.values():
                if self._storage_name.product_id is not None and self._storage_name.product_id != plane.product_id:
                    continue

                for artifact in plane.artifacts.values():
                    if self._storage_name.file_uri != artifact.uri:
                        continue
                    if GemName.is_preview(artifact.uri):
                        continue
                    if self._lookup._reader:
                        update_artifact_meta(artifact, self._lookup._reader.file_info.get(artifact.uri))
                    else:
                        update_artifact_meta(artifact, self._storage_name.file_info.get(artifact.uri))
                    self._update_artifact(artifact)
                    processed = ofr.is_processed(self._storage_name.file_name)
                    if self._instrument in [
                        Inst.MICHELLE,
                        Inst.TRECS,
                        Inst.GNIRS,
                    ]:
                        # Michelle is a retired visitor instrument.
                        # Spatial WCS info is in primary header. There
                        # are a variable number of FITS extensions
                        # defined by primary keyword NUMEXT; assume the
                        # same spatial WCS for each for now - it differs
                        # only slightly because of telescope 'chopping'
                        # and 'nodding' used in acquisition. DB - 01-18-19
                        #
                        # DB - 01-18-19 - GNIRS has no WCS info in extension;
                        # use primary header
                        self._update_position_from_zeroth_header(artifact)

                    delete_these_parts = []
                    for part in artifact.parts:

                        if part == '2' and self._instrument is Inst.GPI:
                            # GPI data sets have two extensions. First is
                            # science image (with WCS), second is data quality
                            # for each pixel (no WCS).
                            self._logger.info(
                                f'GPI: Setting chunks to None for part {part} for {self._observation.observation_id}'
                            )
                            artifact.parts[part].chunks = TypedList(
                                Chunk,
                            )
                            continue
                        for c in artifact.parts[part].chunks:
                            # example is CIRPASS/2003jun30_3385.fits - older
                            # versions of the file had more headers
                            if int(part) >= len(self._headers):
                                delete_these_parts.append(part)
                                continue

                            # energy WCS
                            filter_name = self._get_filter_name(int(part))
                            if self._reset_energy(self._observation.type, filter_name):
                                cc.reset_energy(c)
                            else:
                                self._update_energy(c, plane.data_product_type, filter_name)
                            # position WCS
                            if self._reset_position(self._observation.type, artifact.product_type):
                                self._logger.warning(f'Resetting position for part {part}')
                                cc.reset_position(c)
                            else:
                                self._update_position(part, c, int(part))

                            # time WCS
                            self._update_time(c)

                            self._make_axes_consistent(c)

                    if isinstance(self._observation, DerivedObservation):
                        values = cc.find_keywords_in_headers(self._headers[1:], ['IMCMB'])
                        repaired_values = _remove_processing_detritus(values, self._observation.observation_id)
                        cc.update_plane_provenance_from_values(
                            plane,
                            self._repair_provenance_value,
                            repaired_values,
                            self._storage_name.collection,
                            self._observation.observation_id,
                        )

                    if (
                        processed
                        or isinstance(self._observation, DerivedObservation)
                        or self._instrument is Inst.TEXES
                    ) and 'jpg' not in self._storage_name.file_name:
                        # not the preview artifact
                        if plane.provenance is not None:
                            if self._instrument is not Inst.GRACES:
                                plane.provenance.reference = (
                                    f'http://archive.gemini.edu/searchform/filepre={self._storage_name.file_name}'
                                )

                    for part in delete_these_parts:
                        self._logger.warning(f'Delete part {part} from artifact {artifact.uri}')
                        artifact.parts.pop(part)

                if self._observation.proposal is not None and self._observation.proposal.pi_name is None:
                    program = self._storage_name.md_context.pi_metadata.get(self._observation.proposal.id)
                    if program:
                        self._observation.proposal.pi_name = program.get('pi_name')
                        self._observation.proposal.title = program.get('title')

            if isinstance(self._observation, DerivedObservation):
                cc.update_observation_members(self._observation)

            GeminiMapping.value_repair.repair(self._observation)
        except Exception as e:
            self._logger.error(f'Error {e} for {self._observation.observation_id}')
            tb = traceback.format_exc()
            self._logger.debug(tb)
            raise mc.CadcException(e)
        self._logger.debug('Done update.')
        return self._observation

    def _build_chunk_energy(self, chunk, filter_name):
        if self.fm is not None:
            # If n_axis=1 (as I guess it will be for all but processes GRACES
            # spectra now?) that means crpix=0.5 and the corresponding crval would
            # central_wl - bandpass/2.0 (i.e. the minimum wavelength).   It is fine
            # if you instead change crpix to 1.0.   I guess since the ‘range’ of
            # one pixel is 0.5 to 1.5.

            axis = CoordAxis1D(axis=Axis(ctype='WAVE', cunit='um'))
            ref_coord1 = RefCoord(0.5, self.fm.central_wl - self.fm.bandpass / 2.0)
            ref_coord2 = RefCoord(1.5, self.fm.central_wl + self.fm.bandpass / 2.0)
            axis.range = CoordRange1D(ref_coord1, ref_coord2)

            # DB - 14-02-19 value looks clearer (as two filters) with a space on
            # both sides of the ‘+’.
            bandpass_name = None if (filter_name is None or len(filter_name) == 0) else filter_name.replace('+', ' + ')

            if bandpass_name is not None and 'empty' in bandpass_name:
                # DB 02-12-19 - But some files have filters open1-6 and empty_01.
                # So likely less confusing to remove ‘empty*’ completely.
                #
                # Remove the 'empty' string here, now that min/max wavelength
                # calculations have been completed.
                bandpass_name = bandpass_name.replace('empty', '')
                if len(bandpass_name) == 0:
                    bandpass_name = None
                elif '+' in bandpass_name:
                    bandpass_name = bandpass_name.replace('+', '').strip()

            if math.isclose(self.fm.central_wl, 0.0):
                energy = SpectralWCS(
                    axis=CoordAxis1D(axis=Axis(ctype='WAVE', cunit='um')),
                    specsys='TOPOCENT',
                    ssyssrc=None,
                    ssysobs=None,
                    bandpass_name=bandpass_name,
                )
            else:
                energy = SpectralWCS(
                    axis=axis,
                    specsys='TOPOCENT',
                    ssyssrc='TOPOCENT',
                    ssysobs='TOPOCENT',
                    bandpass_name=bandpass_name,
                    resolving_power=self.fm.resolving_power,
                )
            chunk.energy = energy
            # no chunk energy is derived from FITS file axis metadata, so no
            # cutouts to support
            chunk.energy_axis = None
        else:
            # the goal is findable data, so log the missing metadata, but let the ingestion continue
            self._logger.warning(f'No filter metadata for filter {filter_name}')

    def _get_data_label(self):
        return self._lookup.data_label(self._storage_name.file_uri)

    def _get_filter_name(self, ext):
        """
        Create the filter names for use by update_energy methods.

        :return: The filter names, or None if none found.
        """
        filter_name = self._lookup.filter_name(self._storage_name.file_uri)
        filter_name = self._scrub_filter_name(filter_name)
        if filter_name is None or len(filter_name.strip()) == 0:
            result = self._search_through_keys(ext, ['FILTER'])
            filter_name = result
        self._logger.debug(f'Filter names are "{filter_name}" in {self._storage_name.obs_id}')
        return filter_name

    def _make_axes_consistent(self, chunk):
        # DB 04-17-21
        # BIASes and DARKs for all instruments should ignore spatial and
        # spectral wcs - clean up associated axes
        #
        # also fix a very specific edge case where cal files have useless WCS
        # information for the purposes of CAOM2.4 axis checks, and the
        # corresponding cutouts. No spatial wcs means invalid chunk.naxis
        # value, so set that to None, which then invalidates the
        # chunk.time_axis value DB 06-01-21 - no calibration file cutouts to
        # support, so removing this axis information is not removing
        # downstream functionality
        if chunk.naxis == 3 and chunk.position is None and chunk.time is not None:
            if chunk.time.axis.function.naxis == 1:
                chunk.naxis = None
                chunk.time_axis = None
            else:
                chunk.naxis = 1
                chunk.time_axis = 1

        if (
            (chunk.naxis is not None and chunk.naxis <= 2)
            and not
            # the following exempts the Fox use case
            (chunk.naxis == 1 and chunk.time_axis == 1)
        ):
            if chunk.position_axis_1 is None:
                chunk.naxis = None
            chunk.time_axis = None

    def _multiple_filter_lookup(self, filter_name, lookup, wl_max=None):
        w_max = 10.0 if wl_max is None else wl_max
        w_min = 0.0
        for ii in filter_name.split('+'):
            if ii in lookup:
                wl_max = lookup[ii][2]
                wl_min = lookup[ii][1]
            else:
                msg = f'{self._instrument} Unprepared for filter {ii} from {self._storage_name.obs_id}'
                if self._instrument is Inst.MICHELLE and (ii.startswith('I') or (ii == 'Grid_T')):
                    self._logger.info(msg)
                    continue
                else:
                    raise mc.CadcException(msg)
            if wl_max < w_max:
                w_max = wl_max
            if wl_min > w_min:
                w_min = wl_min
        return w_max, w_min

    def _repair_provenance_value(self, value, obs_id):
        self._logger.debug(f'Being _repair_provenance_value of {value} for {obs_id}.')
        prov_file_id = value
        try:
            uri = mc.build_uri(self._storage_name.collection, prov_file_id, self._storage_name.scheme)
            # TODO - what to do about the data label search here? there's no
            # file metadata, it would be _most_ efficient to do this from CADC
            # so need to re-use this code somehow :(
            prov_obs_id = self._provenance_finder.get(uri)
        except mc.CadcException as e:
            # the file id probably does not exist at on disk, at CADC, or at
            # Gemini, ignore, because it's provenance
            self._logger.warning(f'Failed to find {prov_file_id}')
            # DB 01-06-21 - use not found for the DATALAB/observationID value
            # so it's easy to find in the database and let Gemini know.
            prov_obs_id = 'not_found'
        self._logger.debug(f'End _repair_provenance_value. {prov_obs_id} {prov_file_id}')
        return prov_obs_id, prov_file_id

    def _reset_energy(self, observation_type, filter_name):
        result = False
        om_filter_name = self._lookup.filter_name(self._storage_name.file_uri)
        if (
            observation_type in ['BIAS', 'DARK']
            or (
                self._instrument
                in [
                    Inst.GMOS,
                    Inst.GMOSN,
                    Inst.GMOSS,
                ]
                and observation_type in ['BIAS', 'MASK']
            )
            or (om_filter_name is not None and ('blank' in om_filter_name or 'Blank' in om_filter_name))
            or (filter_name is not None and ('unknown' in filter_name or filter_name == ''))
        ):
            self._logger.info(
                f'No chunk energy for {self._storage_name.obs_id} obs type {observation_type} JSON filter name '
                f'{om_filter_name} FITS filter name {filter_name}.'
            )
            # 'unknown' in filter_name test obs is GN-2004B-Q-30-15-002
            # DB 23-04-19 - GN-2004B-Q-30-15-002: no energy
            # GMOS GS-2005A-Q-26-12-001.  Lots of missing metadata, including
            # release date so no energy (filter_name == '')
            # DB 04-17-21
            # BIASes and DARKs for all instruments should ignore spectral wcs
            result = True
        return result

    def _reset_position(self, observation_type, artifact_type):
        """
        Return True if there should be no spatial WCS information created at
        the chunk level.
        """
        result = False
        types = self._lookup.types(self._storage_name.file_uri)
        ra = self.get_ra(0)
        if ('AZEL_TARGET' in types and ra is None) or observation_type in [
            'BIAS',
            'DARK',
        ]:
            # DB - 02-04-19 - Az-El coordinate frame likely means the
            # telescope was parked or at least not tracking so spatial
            # information is irrelevant.

            # DB - 09-04-19 - AZEL_TARGET should likely be checked for all
            # datasets, and means the spatial WCS should be ignored. since
            # this generally means the telescope is not tracking and so
            # spatial WCS info isn’t relevant since the position is changing
            # with time.

            # DB 04-17-21
            # BIASes and DARKs for all instruments should ignore spatial wcs

            result = True

        # DB 15-03-22
        # ignoring spatial coordinates whenever TRKFRAME or FRAME had the
        # value AZEL_TOPO because the telescope is parked at the zenith.
        trkframe = self._headers[0].get('TRKFRAME')
        frame = self._headers[0].get('FRAME')
        if trkframe is None and frame is None and len(self._headers) > 1:
            # just in case the 0th header is not the source of the metadata
            trkframe = self._headers[1].get('TRKFRAME')
            frame = self._headers[1].get('FRAME')

        if (trkframe is not None and trkframe == 'AZEL_TOPO') or (frame is not None and frame == 'AZEL_TOPO'):
            result = True

        # DB 15-07-22
        # if FRAME is 'No Value' the telescope is not tracking, ignore
        # Spatial WCS for calibration artifacts
        if artifact_type == ProductType.CALIBRATION and frame == 'No Value':
            result = True

        return result

    def _scrub_filter_name(self, filter_name):
        # DB 24-04-19
        # ND = neutral density and so any ND* filter can be ignored as it
        # shouldn’t affect transmission band.  Likely observing a bright
        # Other instruments occasionally have ND filters in the beam.
        if filter_name is not None:
            filter_name = filter_name.replace('&', '+')
            temp = filter_name.split('+')
            for fn in temp:
                if fn.startswith('ND'):
                    filter_name = filter_name.replace(fn, '')
            filter_name = filter_name.strip('+')
        return filter_name

    def _search_through_keys(self, ext, search_keys):
        result = []
        for key in self._headers[ext].keys():
            for lookup in search_keys:
                if lookup in key:
                    value = self._headers[ext].get(key).lower()
                    ignore = False
                    for ii in FILTER_VALUES_TO_IGNORE:
                        if ii.startswith(value) or value.startswith(ii):
                            ignore = True
                            break
                    if ignore:
                        continue
                    else:
                        result.append(self._headers[ext].get(key).strip())
        return '+'.join(result)

    def _update_artifact(self, artifact):
        pass

    def _update_composite(self, obs):
        result = None
        if self._instrument is Inst.TRECS:
            if self._storage_name.product_id is not None and (
                self._storage_name.product_id.startswith('rS') or self._storage_name.product_id.startswith('rN')
            ):
                # DB 02-06-20
                # processed TReCS files in Gemini's archive are derived by
                # combining the NNODSETS x NSAVSETS contained within a single
                # unprocessed image into a simpler image array.
                #
                # SGo - this means ignoring the IMCMB keywords that are an
                # artifact of that, which is how Composite construction is
                # otherwise determined.
                result = obs
        else:
            result = cc.change_to_composite(obs)
            self._logger.info(f'{obs.observation_id} is a Derived Observation.')
        return result

    def _update_position(self, part, chunk, extension):
        self._logger.debug(f'Enter _update_position for {self._storage_name.file_uri}')
        if chunk.position is not None:
            # DB 14-07-22
            # test if the absolute value of any CD element is > 1/3600
            # (one arcsecond/pixel). This is much larger than any Gemini
            # pixel would be and so the spatial WCS would likely be nonsense
            # for any science OR calibration image.
            fail = False
            for keyword in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                temp = mc.get_keyword(self._headers, keyword)
                if temp is not None:
                    if abs(temp) > 1.0 / 3600.0:
                        fail = True
                        break

            if fail:
                self._logger.warning(f'Large CD matrix values, no Spatial WCS for {self._storage_name.file_uri}')
                cc.reset_position(chunk)
        self._logger.debug('End _update_position')

    def _update_position_from_zeroth_header(self, artifact):
        """Make the 0th header spatial WCS the WCS for all the
        chunks."""
        primary_header = self._headers[0]

        # naxis values are only available from extensions

        primary_header['NAXIS1'] = self._headers[-1].get('NAXIS1')
        primary_header['NAXIS2'] = self._headers[-1].get('NAXIS2')

        primary_chunk = Chunk()
        types = self._lookup.types(self._storage_name.file_uri)
        ra = self._headers[0].get('RA')
        ra_json = self.get_ra(0)
        dec_json = self._lookup.dec(self._storage_name.file_uri)
        if ('AZEL_TARGET' in types and ra is None) or (
            self._instrument is Inst.GNIRS and ra_json is None and dec_json is None
        ):
            # DB - 02-04-19 - Az-El coordinate frame likely means the telescope
            # was parked or at least not tracking so spatial information is
            # irrelevant.

            # DB - 09-04-19 - AZEL_TARGET should likely be checked for all
            # datasets, and means the spatial WCS should be ignored. since this
            # generally means the telescope is not tracking and so spatial WCS
            # info isn’t relevant since the position is changing with time.

            # DB 24-04-19
            # Ignore spatial WCS if RA/Dec are not available for GNIRS.

            self._logger.info(f'Spatial WCS is None for {self._storage_name.obs_id}')
        else:
            wcs_parser = FitsWcsParser(primary_header, self._storage_name.obs_id, 0)
            wcs_parser.augment_position(primary_chunk)

        for part in artifact.parts:
            if part == '0':
                continue
            for chunk in artifact.parts[part].chunks:
                if primary_chunk.position is not None:
                    chunk.position = primary_chunk.position
                    chunk.position_axis_1 = 1
                    chunk.position_axis_2 = 2

    def _update_time(self, chunk):
        # the default is to do nothing
        pass


def _remove_processing_detritus(values, obs_id):
    """
    There are several naming patterns in the provenance for processed files.
    Try to extract meaningful raw file names.

    Remove duplicates, so values used for provenance in multiple forms are
    only looked up once.
    :param values: A list of IMCMB* keyword values.
    :param obs_id: str for logging information
    :return: A list of unique provenance file names, such as may be found
        at Gemini. The list may be empty.
    """
    logging.debug(f'Begin _remove_processing_detritus for {obs_id}')
    result = []
    for value in values:
        # e.g.
        # IMCMB001 = 'tmpimgwsk9476kd_5.fits[SCI,1]'
        # tmpfile22889S20141226S0203.fits[SCI,1]
        # IMCMB001= 'rawdir$2004may20_0048.fits'

        if 'N' in value:
            temp = 'N' + value.split('N', 1)[1]
        elif 'S' in value:
            x = value.split('S')
            if len(x) == 2 and '[SCI' in value:
                logging.warning(f'Unrecognized IMCMB value {value}')
                continue
            else:
                temp = 'S' + value.split('S', 1)[1]
        elif '$' in value:
            temp = value.split('$', 1)[1]
        else:
            logging.warning(f'Unrecognized IMCMB value {value}')
            continue

        if '_' in temp and (temp.startswith('S') or temp.startswith('N')):
            temp1 = temp.split('_')[0]
        elif '.fits' in temp:
            temp1 = temp.split('.fits')[0]
        elif '[SCI' in temp:
            temp1 = temp.split('[SCI')[0]
        else:
            logging.warning(f'Failure to repair {temp}')
            continue

        result.append(temp1[:14])

    logging.debug('End _remove_processing_detritus')
    return list(set(result))


class Bhros(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')

    def get_dec(self, ext):
        # bHROS, TEXES ra/dec not in json
        return self._headers[ext].get('DEC')

    def get_ra(self, ext):
        # bHROS, TEXES: ra/dec not in json
        return self._headers[ext].get('RA')

    def _reset_energy(self, observation_type, filter_name):
        return observation_type in ['BIAS', 'DARK']

    def _update_energy(self, chunk, data_product_type, filter_name):
        """bhros-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy')
        # DB - 20-02-19 - There were bHROS filters but I don’t think they were
        # used during the very limited lifetime of the instrument.  No info
        # in the headers either.
        #
        # bHROS spectral resolution should be approximately 150,000/x-binning
        # value. json returns a “detector_binning”: “1x1" value where the 1x1
        # indicates no binning in x or y (for this example).  Could be 2x1,
        # 4x2, etc.  Binning is determined from header keyword:
        # CCDSUM  = ‘1 1     ’           / CCD pixel summing
        #
        # The approximate central wavelength is json ‘central_wavelength’
        # value. Unfortunately the wavelength coverage is not straightforward.
        # See
        # http://www.gemini.edu/sciops/instruments/hros/hrosDispersion.html.
        # The CCD did not cover the entire spectrum so only a subset of the
        # entire optical spectral region was observed.

        # Use central wavelength in microns and +/- 0.2 microns as a better
        # guess-timate rather than imply that entire spectrum is present.

        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS spectroscopy for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = self._lookup.central_wavelength(self._storage_name.file_uri)
            self.fm.adjust_bandpass(0.2)
            self.fm.resolving_power = 150000.0
            ccd_sum = self._headers[0].get('CCDSUM')
            if ccd_sum is not None:
                temp = float(ccd_sum.split()[1])
                self.fm.resolving_power = 150000.0 / temp
            self._build_chunk_energy(chunk, filter_name)
        else:
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._logger.debug(f'End _update_energy')

    def _update_position(self, part, chunk, ext):
        self._headers[0]['CTYPE1'] = 'RA---TAN'
        self._headers[0]['CTYPE2'] = 'DEC--TAN'
        self._headers[0]['CUNIT1'] = 'deg'
        self._headers[0]['CUNIT2'] = 'deg'
        self._headers[0]['CRVAL1'] = self.get_ra(0)
        self._headers[0]['CRVAL2'] = self.get_dec(0)
        self._headers[0]['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        self._headers[0]['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        self._headers[0]['CROTA1'] = 0.0
        self._headers[0]['NAXIS1'] = 1
        self._headers[0]['NAXIS2'] = 1
        self._headers[0]['CRPIX1'] = self.get_crpix1(0)
        self._headers[0]['CRPIX2'] = self.get_crpix2(0)
        self._headers[0]['CD1_1'] = self.get_cd11(0)
        self._headers[0]['CD1_2'] = 0.0
        self._headers[0]['CD2_1'] = 0.0
        self._headers[0]['CD2_2'] = self.get_cd22(0)
        self._headers[0]['EQUINOX'] = mc.to_float(self._headers[0].get('TRKEQUIN'))
        wcs_parser = FitsWcsParser(self._headers[0], self._storage_name.obs_id, 0)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        chunk.position.coordsys = self._headers[0].get('TRKFRAME')
        self._logger.debug('End _update_chunk_position')


class Cirpass(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')
        bp.configure_position_axes((1, 2))

    def get_cd11(self, ext, keyword='CDELT1'):
        lens_scl = self._headers[ext].get('LENS_SCL')
        return mc.to_float(lens_scl) / 3600.0

    def get_cd22(self, ext):
        return self.get_cd11(ext)

    def get_data_product_type(self, ext):
        return DataProductType.SPECTRUM

    def get_dec(self, ext):
        # DB - 06-03-19 - Must use FITS header info for most WCS info
        ra, dec = self.get_sky_coord(ext, 'TEL_RA', 'TEL_DEC')
        return dec

    def get_exposure(self, ext):
        # exposure_time is null in the JSON
        # DB - 06-03-19 Use value of header keyword EXP_TIME for exposure time.
        return mc.to_float(self._headers[ext].get('EXP_TIME'))

    def get_obs_intent(self, ext):
        data_label = self._headers[ext].get('DATALAB')
        if data_label is not None and '-CAL' in data_label:
            result = ObservationIntentType.CALIBRATION
        else:
            obs_type = self._headers[ext].get('OBS_TYPE')
            if obs_type in CAL_VALUES:
                result = ObservationIntentType.CALIBRATION
            else:
                result = ObservationIntentType.SCIENCE
        return result

    def get_obs_type(self, ext):
        # DB - 06-03-19
        temp = self._headers[ext].get('OBJECT')
        # if ‘dome’ or ‘flat’ or ‘twilight’ or ‘gcal’ in
        #       object string:  obstype = FLAT
        # if ‘argon’ or ‘arc’ in string: obstype = ARC
        # if ‘dark’ in string: obstype = DARK
        # Otherwise assume obstype = OBJECT
        result = 'OBJECT'
        for ii in ['dome', 'flat', 'twilight', 'gcal']:
            if ii in temp:
                result = 'FLAT'
                break
        for ii in ['argon', 'arc']:
            if ii in temp:
                result = 'ARC'
                break
        if 'dark' in temp:
            result = 'DARK'
        return result

    def get_ra(self, ext):
        # DB - 06-03-19 - Must use FITS header info for most WCS info
        ra, dec = self.get_sky_coord(ext, 'TEL_RA', 'TEL_DEC')
        return ra

    def _get_filter_name(self, ext):
        # filter names are only in HDU0?
        return super()._get_filter_name(0)

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        filter_name = ''
        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS spectral mode for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata()
            self.fm.set_central_wl(1.0, 1.67)
            self.fm.set_bandpass(1.0, 1.67)
            self.fm.resolving_power = 3200.0
        else:
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_position(self, part, chunk, ext):
        self._headers[0]['CTYPE1'] = 'RA---TAN'
        self._headers[0]['CTYPE2'] = 'DEC--TAN'
        self._headers[0]['CUNIT1'] = 'deg'
        self._headers[0]['CUNIT2'] = 'deg'
        self._headers[0]['CRVAL1'] = self.get_ra(0)
        self._headers[0]['CRVAL2'] = self.get_dec(0)
        self._headers[0]['CROTA1'] = 0.0
        # So perhaps try:
        #     NAXIS1 = 33
        #     NAXIS2 = 15
        self._headers[0]['NAXIS1'] = 33
        self._headers[0]['NAXIS2'] = 15

        # TODO TODO TODO - add the get_crpix1, get_crpix2, get_cd11, get_cd22
        # methods, because that may allow for extraction of common code??????

        # LENS_SCL determines the scale/lenslet:  0.36 or 0.25 (arcseconds
        # per lens)
        #     if 0.36 then FOV is 13.0" x 4.7" (RA and Dec)
        #     if 0.25 then FOV is 9.3" x 3.5"
        lens_scl = self._headers[0].get('LENS_SCL')
        if lens_scl == '0.36':
            self._headers[0]['CDELT1'] = 13.0 / 3600.0
            self._headers[0]['CDELT2'] = 4.7 / 3600.0
        else:
            self._headers[0]['CDELT1'] = 9.3 / 3600.0
            self._headers[0]['CDELT2'] = 3.5 / 3600.0
        # DB - 06-03-19
        # FOV is fixed at two possible values and has no bearing on NAXIS1/2
        # values. See
        # http://www.gemini.edu/sciops/instruments/cirpass/cirpassIFU.html.
        # 499 lenslets cover the FOV: about 33 along one axis and 15 along the
        # other.
        #
        # LENS_SCL determines the scale/lenslet:  0.36 or 0.25
        # (arcseconds per lens)
        #     if 0.36 then FOV is 13.0" x 4.7" (RA and Dec)
        #     if 0.25 then FOV is 9.3" x 3.5"
        #     cd11 = LENS_SCL/3600.0
        #     cd22 = LENS_SCL/3600.0
        self._headers[0]['CRPIX1'] = self.get_crpix1(0)
        self._headers[0]['CRPIX2'] = self.get_crpix2(0)
        self._headers[0]['CD1_1'] = self.get_cd11(0)
        self._headers[0]['CD1_2'] = 0.0
        self._headers[0]['CD2_1'] = 0.0
        self._headers[0]['CD2_2'] = self.get_cd22(0)
        wcs_parser = FitsWcsParser(self._headers[0], self._storage_name.obs_id, 0)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        self._logger.debug('End _update_position')


class F2(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        if self.get_obs_type(0) != 'FLAT' or (
            self.get_data_product_type(0) == DataProductType.SPECTRUM and self.get_obs_type(0) == 'FLAT'
        ):
            # WF 07-12-22
            # - be specific about the problem and ingest as best as possible
            # - would like to know when the problem gets fixed (so, need to know when the problem first occurs)
            ctype1 = mc.get_keyword(self._headers, 'CTYPE1')
            ctype2 = mc.get_keyword(self._headers, 'CTYPE2')
            if ctype1 is None or ctype2 is None or '?' in ctype1 or '?' in ctype2:
                self._observable.rejected.record(mc.Rejected.BAD_METADATA, self._storage_name.file_name)
            else:
                bp.configure_position_axes((1, 2))

    def get_obs_type(self, ext):
        # DB 03-06-21
        # check the 'types' JSON value
        result = super().get_obs_type(ext)
        types = self._lookup.types(self._storage_name.file_uri)
        if 'DARK' in types:
            result = 'DARK'
        return result

    def _get_filter_name(self, ext):
        """
        Look up and manipulate into an expected pattern the filter names for
        use by _update_energy methods.

        :return: The filter names, or None if none found.
        """
        for filter_name in [
            self._lookup.filter_name(self._storage_name.file_uri),
            self._search_through_keys(ext, ['FILTER']),
            # OO/WF - 12-07-22
            # Oliver stated that the header keyword “LYOT” is also used to
            # confirm the filter keyword. In the case of S20220203S0019,
            # LYOT=‘J-lo_G0801’ which is a valid filter element.
            #
            # In the cases of FILTER='OPEN' FILTER1='OPEN' FILTER2='OPEN', the
            # logic is to look at LYOT
            self._headers[0].get('LYOT'),
        ]:
            filter_name = self._scrub_filter_name(filter_name)
            if filter_name is None or len(filter_name.strip()) == 0 or filter_name.lower() == 'open':
                continue
            else:
                break

        self._logger.debug(f'Filter names are "{filter_name}" in {self._storage_name.obs_id}')
        return filter_name

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        # DB - 02-05-19
        # For F2 use SVO filter service for bandpass or ‘delta’ for images and
        # spectroscopy.  Treat images as for other instruments but…
        #
        # if gemini_md[‘spectroscopy’] == ‘true’:
        #  ref_wl = gemini_md[‘central_wavelength’]   or  GRWLEN header value
        #  grism = gemini_md[‘disperser’]  or  GRISM header value
        #  if gemini_md[‘mode’] == ‘LS’:  # long-slit
        #    slit_width = MASKNAME header value, e.g. ‘4pix-slit’, but need
        #                                                          only ‘4’
        #    use the table I sent you with slit/grism values to determine
        #    average resolution R
        #
        #  elif gemini_md[‘mode’] == ‘MOS’:  # Multi-object
        #    slit_width = 2   # no way to determine slit widths used in
        #    custom mask, so assume 2
        #    use the table I sent you with slit/grism values to determine
        #    average resolution R
        #
        #  else:
        #    fail because there shouldn’t be any other spectroscopy mode

        # DB 09-04-19 - Ignore energy when the grism is in the header but
        # object value of “COVER CLOSED” so is another type of calibration
        # exposure apparently.

        reset_energy = False
        object_value = self._lookup.object(self._storage_name.file_uri)
        if 'COVER CLOSED' in object_value or 'Undefined' in filter_name:
            # DB 30-04-19
            # Flamingos ‘Undefined’ filter:  no spectral WCS
            reset_energy = True
        else:
            filter_md = self._filter_cache.get_filter_metadata(
                self._instrument,
                filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )
            if data_product_type == DataProductType.IMAGE:
                self._logger.debug(f'SpectralWCS: imaging mode for {self._storage_name.obs_id}.')
                if filter_md is None:
                    self._logger.warning(
                        f'No filter metadata found for {filter_name} in {self._storage_name.file_uri}'
                    )
                    reset_energy = True
                self.fm = filter_md
            elif data_product_type == DataProductType.SPECTRUM:
                self._logger.debug(f'SpectralWCS: LS|Spectroscopy mode for {self._storage_name.obs_id}.')
                fp_mask = self._headers[0].get('MASKNAME')
                mode = self._lookup.mode(self._storage_name.file_uri)
                slit_width = None
                if mode == 'LS':
                    slit_width = fp_mask[0]
                self.fm = svofps.FilterMetadata()
                self.fm.central_wl = filter_md.central_wl
                self.fm.bandpass = filter_md.bandpass
                grism_name = self._headers[0].get('GRISM')
                self._logger.debug(f'grism name is {grism_name} fp_mask is {fp_mask} for {self._storage_name.obs_id}')
                # lookup values from
                # https://www.gemini.edu/sciops/instruments/flamingos2/spectroscopy/longslit-spectroscopy
                lookup = {
                    '1': [1300.0, 3600.0],
                    '2': [900.0, 2800.0],
                    '3': [600.0, 1600.0],
                    '4': [350.0, 1300.0],
                    '6': [130.0, 1000.0],
                    '8': [100.0, 750.0],
                }
                if slit_width is None or slit_width not in lookup:
                    # DB 02-04-19
                    # For F2 at line 1409 of main_app.py set slit_width = ‘2’
                    # as a default of slit_width[0] is not a numeric value
                    slit_width = '2'
                if grism_name.startswith('R3K_'):
                    self.fm.resolving_power = lookup[slit_width][1]
                else:
                    self.fm.resolving_power = lookup[slit_width][0]
            else:
                raise mc.CadcException(
                    f'{self._instrument}: Do not understand DataProductType '
                    f'{data_product_type} for {self._storage_name.obs_id}'
                )

        if reset_energy:
            self._logger.info(f'Setting spectral WCs to none for {self._storage_name.obs_id}')
            cc.reset_energy(chunk)
        else:
            temp = (
                filter_name
                if (filter_md is None or len(filter_md.filter_name) == 0)
                else '+'.join(ii for ii in filter_md.filter_name)
            )
            self._build_chunk_energy(chunk, temp)
        self._logger.debug('End _update_energy')

    def _update_time(self, chunk):
        """F2 FITS files have a CD3_3 element that's not supported by
        fits2caom2, so using the blueprint will not work to adjust that
        value. Set delta specifically here."""
        self._logger.debug(f'Begin _update_time')
        if chunk.time is not None and chunk.time.axis is not None and chunk.time.axis.function is not None:
            exposure = mc.to_float(self._lookup.exposure_time(self._storage_name.file_uri))
            chunk.time.axis.function.delta = mc.convert_to_days(exposure)
            self._logger.info(f'Updated time delta for {self._storage_name.obs_id} to {exposure}.')
        self._logger.debug(f'End _update_time')


class Flamingos(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))
        # DB 27-05-19
        # Flamingos, you actually want to use the EQUINOX value, not the
        # EPOCH.   And I think EQUINOX header value is usually 2000.0, even
        # for the example GS-CAL20020620-15-0462 02jun20.0462 with
        # RA_TEL = “UNAVAILABLE”.  For Gemini the assumption is that the
        # RA/Dec values in the headers are always based on the position of
        # the equinox given at the time specified by the EQUINOX keyword
        # value.
        bp.clear('Chunk.position.equinox')
        bp.add_attribute('Chunk.position.equinox', 'EQUINOX')

    def get_art_product_type(self, ext):
        # DB - 28-02-19
        # Another relatively minor thing for FLAMINGOS:  get_obs_intent likely
        # has to have FLAMINGOS added to it.  If obs_type is DARK,
        # ACQUISIITON, FLAT or ARC then the intent is calibration, otherwise
        # science. cal_values should have those obs_type values in it as well.
        # So if ‘CAL’ is not in the data_label and obs_type is not in
        # cal_values then it’s science.
        object_value = self._lookup.object(self._storage_name.file_uri)
        for ii in ['Dark', 'DARK', 'Arc', 'flat', 'Flat', 'Acq', 'acq']:
            if ii in object_value:
                return ProductType.CALIBRATION
        data_label = self._get_data_label()
        if data_label is not None and 'CAL' in data_label:
            return ProductType.CALIBRATION
        return ProductType.SCIENCE

    def get_data_product_type(self, ext):
        # DB - 18-02-19
        # Determine mode from DECKER and GRISM keywords:
        #
        # If DECKER = imaging
        # observing mode = imaging
        # else if DECKER = mos or slit
        # if GRISM contains the string ‘open’
        # observation mode = imaging and observation type = ACQUISITION
        # else if GRISM contains the string ‘dark’
        # observation mode = spectroscopy and observation type = DARK
        # else
        # observation mode = spectroscopy
        # else
        # error? DECKER might also be ‘dark’ at times in which case
        # observation type = dark and observation mode can be either
        # spectroscopy or imaging but I haven’t found an example yet.
        #
        # As indicated above in if/else code, if GRISM  keyword
        # contains substring ‘dark’ then that is another way to
        # identify a dark observation.
        decker = self._headers[ext].get('DECKER')
        data_type = DataProductType.SPECTRUM
        if decker is None:
            raise mc.CadcException(
                f'{self._instrument}: No mode information found for {self._storage_name.file_name}'
            )
        else:
            if decker == 'imaging':
                data_type = DataProductType.IMAGE
            elif decker in ['mos', 'slit']:
                grism = self._headers[ext].get('GRISM')
                if grism is not None:
                    if 'open' in grism:
                        data_type = DataProductType.IMAGE
                    elif 'dark' in grism:
                        data_type = DataProductType.SPECTRUM
        return data_type

    def get_exposure(self, ext):
        # exposure_time is null in the JSON
        # DB - 06-03-19 Use value of header keyword EXP_TIME for exposure
        # time.
        return mc.to_float(self._headers[ext].get('EXP_TIME'))

    def get_obs_type(self, ext):
        # DB - 18-02-19
        # Determine mode from DECKER and GRISM keywords:
        #
        # If DECKER = imaging
        # observing mode = imaging
        # else if DECKER = mos or slit
        # if GRISM contains the string ‘open’
        # observation mode = imaging and observation type = ACQUISITION
        # else if GRISM contains the string ‘dark’
        # observation mode = spectroscopy and observation type = DARK
        # else
        # observation mode = spectroscopy
        # else
        # error? DECKER might also be ‘dark’ at times in which case
        # observation type = dark and observation mode can be either
        # spectroscopy or imaging but I haven’t found an example yet.
        #
        # As indicated above in if/else code, if GRISM  keyword
        # contains substring ‘dark’ then that is another way to
        # identify a dark observation.
        decker = self._headers[ext].get('DECKER')
        data_type = DataProductType.SPECTRUM
        obs_type = None
        if decker is not None:
            if decker == 'imaging':
                data_type = DataProductType.IMAGE
            elif decker in ['mos', 'slit']:
                grism = self._headers[ext].get('GRISM')
                if grism is not None:
                    if 'open' in grism:
                        data_type = DataProductType.IMAGE
                        obs_type = 'ACQUISITION'
                    elif 'dark' in grism:
                        data_type = DataProductType.SPECTRUM
                        obs_type = 'DARK'
        else:
            raise mc.CadcException(
                f'{self._instrument}: No mode information found for {self._storage_name.file_name}'
            )
        if obs_type is None:
            # DB - Also, since I’ve found FLAMINGOS spectra if OBJECT keyword
            # or json value contains ‘arc’ (any case) then it’s an ARC
            # observation type
            #
            # For FLAMINGOS since OBS_TYPE seems to always be set to ‘Object’
            # could in principal look at the OBJECT keyword value or json
            # value: if it contains ‘flat’ as a substring (any case) then set
            # observation type to ‘flat’.  Ditto for ‘dark’
            object_value = self._lookup.object(self._storage_name.file_uri)
            if object_value is None:
                object_value = self._headers[ext].get('OBJECT')
            object_value = object_value.lower()
            if 'arc' in object_value:
                obs_type = 'ARC'
            elif 'flat' in object_value:
                obs_type = 'FLAT'
            elif 'dark' in object_value:
                obs_type = 'DARK'
            else:
                # DB - 27-02-19 - Default to 'OBJECT'
                obs_type = 'OBJECT'
        return obs_type

    def get_time_function_val(self, ext):
        # Another FLAMINGOS correction needed:  DATE-OBS in header and
        # json doesn’t include the time  but sets it to 00:00:00.  You have
        # two choices:  concatenate DATE-OBS and UTC header values to the
        # standard form “2002-11-06T07:06:00.3” and use as you use json
        # value for computing temporal WCS data or, for FLAMINGOS only,
        # use the MJD header value in the header as the CRVAL for time.
        return self._headers[ext].get('MJD')

    def _reset_position(self, observation_type, artifact_type):
        result = super()._reset_position(observation_type, artifact_type)
        ra_tel = self._headers[0].get('RA_TEL')
        if ra_tel == 'Unavailable':
            result = True
        return result

    def _update_energy(self, chunk, data_product_type, filter_name):
        """Flamingos-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin _update_energy')

        # DB - 18-02-19 - FLAMINGOS spectral WCS should be similar to what was
        # done for NIRI.  Use NAXIS1 keyword value to determine number of
        # pixels.  The GRISM and FILTER keywords give the same filter ID.
        # The SVO filter information for KPNO/Flamingos has a blue leak in the
        # JK blocking filter which give too large a spectral range.  Can you
        # instead hard code min/max wavelengths using the top two lines in
        # this table 7:
        # http://www-kpno.kpno.noao.edu/manuals/flmn/flmn.user.html#flamspec.
        # Use these for min/max wavelengths and use the average as the
        # ‘central’ wavelength.  Spectral resolution is fixed at about 1300
        # for both grisms.

        # spectral wcs units are microns, values from Table 7 are angstroms.
        # The conversion is here.

        # DB - 21-02-19 - JH central=1.45 um, FWHM=0.95 um
        # 0 = central wavelength
        # 1 = FWHM
        lookup = {
            'JH': [1.45, 0.95],
            'HK': [(2.7588 + 1.0347) / 2.0, (2.7588 - 1.0347)],
        }
        self.fm = svofps.FilterMetadata()
        if filter_name in lookup:
            self.fm.central_wl = lookup[filter_name][0]
            self.fm.bandpass = lookup[filter_name][1]
        else:
            self.fm = self._filter_cache.get_filter_metadata(
                self._instrument,
                filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )
            if self.fm is None:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(
                    f'{self._instrument}: Mystery filter {filter_name} for {self._storage_name.obs_id}'
                )
        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS for {self._storage_name.obs_id}.')
            self.fm.resolving_power = 1300.0
        elif data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging mode for {self._storage_name.obs_id}.')
        else:
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug(f'End _update_energy')

    def _update_position(self, part, chunk, ext):
        # DB - I see nothing in astropy that will do a transformation from
        # crota form to CD matrix, but this is it:

        # cd1_1 = cdelt1 * cos (crota1)
        # cd1_2 = -cdelt2 * sin (crota1)
        # cd2_1 = cdelt1 * sin (crota1)
        # cd2_2 = cdelt2 * cos (crota1)

        # Note that there is not a crota2 keyword (it would have the same
        # value as crota1 if it existed)
        if (
            chunk is not None
            and chunk.position is not None
            and chunk.position.axis is not None
            and chunk.position.axis.function is not None
        ):
            header = self._headers[ext]
            crval1 = header.get('CRVAL1')
            crval2 = header.get('CRVAL2')
            if 0.0 <= crval1 <= 360.0 and -90.0 <= crval2 <= 90.0:
                c_delt1 = header.get('CDELT1')
                c_delt2 = header.get('CDELT2')
                c_rota1 = header.get('CROTA1')
                if c_delt1 is not None and c_delt2 is not None and c_rota1 is not None:
                    chunk.position.axis.function.cd11 = c_delt1 * math.cos(c_rota1)
                    chunk.position.axis.function.cd12 = -c_delt2 * math.sin(c_rota1)
                    chunk.position.axis.function.cd21 = c_delt1 * math.sin(c_rota1)
                    chunk.position.axis.function.cd22 = c_delt2 * math.cos(c_rota1)
                else:
                    self._logger.info(f'Missing spatial wcs inputs for {self._storage_name.obs_id}')
                    chunk.position.axis.function.cd11 = None
                    chunk.position.axis.function.cd12 = None
                    chunk.position.axis.function.cd21 = None
                    chunk.position.axis.function.cd22 = None
            else:
                # DB 04-12-19
                # FLAMINGOS GS-CAL20020623-14-0080 02jun23.0080.fits
                # The header has “CRVAL1  =           3581.13808 ” which is
                # supposed to be in degrees and shouldn’t be > 360. Skip
                # spatial WCS if errors like this occur.
                self._logger.warning(
                    f'Spatial WCS set to None for {self._storage_name.obs_id} '
                    f'because CRVAL1 == {crval1} and CRVAL2 == {crval2}.'
                )
                cc.reset_position(chunk)
        else:
            self._logger.info(f'Missing spatial wcs for {self._storage_name.obs_id}')


class Fox(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.set(
            'Plane.provenance.reference',
            f'http://archive.gemini.edu/searchform/filepre={self._storage_name.file_id}.fits',
        )
        bp.clear('Chunk.time.axis.function.naxis')
        bp.add_attribute('Chunk.time.axis.function.naxis', 'NAXIS3')
        bp.set_default('Chunk.time.axis.function.naxis', 1)
        bp.configure_position_axes((1, 2))
        bp.set('Chunk.position.equinox', 'get_position_equinox()')

    def get_data_product_type(self, ext):
        # DB 31-08-20
        # Both cameras are used for speckle imaging.  Datasets consist of
        # cubes of 256 x 256 x 1000 images.  i.e. 1000 short exposures << 1
        # second long with 256 x 256 pixel images (or smaller images if the
        # CCD is binned).
        #
        # So dataProductType = cube
        #
        # PD 02-09-20
        # if those two files are images in the normal sense then it could make
        # sense to create separate planes with dataProductType = image that
        # end up with the correct (distinct) energy metadata.
        return DataProductType.IMAGE

    def get_position_equinox(self, ext):
        result = self._headers[ext].get('EQUINOX', self._headers[ext].get('EPOCH'))
        if result and 1800.0 <= result <= 2500:
            return result
        else:
            return None

    def get_target_type(self, ext):
        return TargetType.OBJECT

    def _update_energy(self, chunk, data_product_type, filter_name):
        """General chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin _update_energy')
        if data_product_type is DataProductType.IMAGE:
            self._logger.debug(f'Spectral WCS {data_product_type} mode for {self._storage_name.obs_id}.')
            self.fm = self._filter_cache.get_filter_metadata(
                self._instrument,
                filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )
            if self.fm is None:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(f'{self._instrument}: mystery filter {filter_name}')
            self._build_chunk_energy(chunk, filter_name)
        else:
            raise mc.CadcException(
                f'{self._instrument} no Spectral WCS support when '
                f'DataProductType {data_product_type} for '
                f'{self._storage_name.obs_id}'
            )
        self._logger.debug(f'End _update_energy')

    def _update_time(self, chunk):
        """
        DB 02-09-20
        Exposure time using JSON values isn’t correct.  I know that for this
        example Gemini shows the exposure time is 0.02 seconds but there are
        1000 x 0.02-second exposures in the cube.  The keyword EXPOSURE gives
        the total exposure time (in seconds), time.exposure, or 20 in this
        case while the json exposure_time should be the time.resolution.
        """
        self._logger.debug(f'Begin _update_time')
        if chunk.time is not None:
            chunk.time.exposure = self._headers[0].get('EXPOSURE')
            # chunk.time.resolution already set by blueprint
        self._logger.debug(f'End _update_time')


class GHOSTSpectralTemporal(GeminiMapping):
    """Fibre-fed spectrograph"""

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super(GeminiMapping, self).accumulate_blueprint(bp)
        self._accumulate_obs_plane_artifact_blueprint(bp)
        # DB 02-02-24
        # for temporal WCS use JSON ut_datetime for the starting time, exposure_time for the exposure time and add
        # the exposure_time to the ut_datetime to determine the end time for ALL obstypes.
        self._accumulate_chunk_time_axis_blueprint(bp, 4)
        bp.configure_energy_axis(3)
        for extension in range(0, len(self._headers)):
            bp.set('Chunk.energy.specsys', 'TOPOCENT', extension)

            camera = self._headers[extension].get('CAMERA')
            if camera and camera in ['RED', 'BLUE']:
                # if camera is SLITV, there's no interesting science information in the HDU
                bp.set('Chunk.energy.axis.axis.ctype', 'WAVE', extension)
                bp.set('Chunk.energy.axis.axis.cunit', 'nm', extension)
                bp.set('Chunk.energy.axis.range.start.pix', 0.5, extension)
                bp.set('Chunk.energy.axis.range.start.val', '_get_energy_chunk_range_start_val()', extension)
                bp.set('Chunk.energy.axis.range.end.pix', 1.5, extension)
                bp.set('Chunk.energy.axis.range.end.val', '_get_energy_chunk_range_end_val()', extension)
                bp.set('Chunk.energy.resolvingPower', '_get_energy_chunk_resolving_power()', extension)

        self._logger.debug(f'End accumulate_blueprint')

    def get_target_type(self, ext):
        return TargetType.OBJECT

    def _get_energy_chunk_range_end_val(self, ext):
        # DB 18-01-24
        # blue channel wavelength coverage is 347 to 530 nm; red channel is 520 to 1060 nm. So 347 to 1060 nm total.
        result = None
        camera = self._headers[ext].get('CAMERA')
        if camera == 'RED':
            result = 1060.0  # nm
        elif camera == 'BLUE':
            result = 530.0  # nm
        return result

    def _get_energy_chunk_range_start_val(self, ext):
        # DB 18-01-24
        # blue channel wavelength coverage is 347 to 530 nm; red channel is 520 to 1060 nm. So 347 to 1060 nm total.
        result = None
        camera = self._headers[ext].get('CAMERA')
        if camera == 'RED':
            result = 520.0  # nm
        elif camera == 'BLUE':
            result = 347.0  # nm
        return result

    def _get_energy_chunk_resolving_power(self, ext):
        """
        DB 18-01-24
        1. resolution is determined by the value for the PHU RESOLUT keyword.  If the value is Standard then the
        starting point for the resolution is 56,000; if High then it's 76,000.

        2. PHU keywords REDCCDS and BLUCCDS give the binning factor along the X (dispersion) and Y (spatial)
        direction.  e.g. a value of 2 x 4 would reduce the resolution by a factor of about 2 (the first number
        in that value).
        """
        result = None
        start = None
        resolut = self._headers[0].get('RESOLUT')
        if resolut == 'Standard':
            start = 56000.0
        elif resolut == 'High':
            start = 76000.0
        elif resolut == 'No Value':
            return result

        if start:
            red_ccds = self._headers[0].get('REDCCDS')
            blue_ccds = self._headers[0].get('BLUCCDS')

            factor_1 = _to_float(red_ccds.split('x')[0])
            factor_2 = _to_float(blue_ccds.split('x')[0])
            if factor_1 == factor_2:
                factor = factor_1
            else:
                raise mc.CadcException(
                    f'Factors do not match for {self._storage_name.file_uri}. Blue {factor_2}, Red {factor_1}'
                )
            result = start / factor
        else:
            self._logger.error(f'Could not determine resolution for {resolut} for {self._storage_name.file_uri}')
        return result

    def _reset_energy(self, observation_type, filter_name):
        return False

    def _update_energy(self, part, chunk, extension):
        pass

    def _reset_position(self, observation_type, artifact_type):
        return

    def _update_position(self, part, chunk, extension):
        pass

    def _update_artifact(self, artifact):
        self._logger.debug(f'Begin _update_artifact for {artifact.uri}')
        for part in artifact.parts.values():
            extension = int(part.name)
            camera = self._headers[extension].get('CAMERA')
            if camera and camera == 'SLITV':
                part.product_type = ProductType.AUXILIARY
        self._logger.debug('End _update_artifact')

    def _update_time(self, chunk):
        if (
            chunk is not None
            and chunk.time is not None
            and chunk.time.axis is not None
            and chunk.time.axis.function is not None
            and chunk.time.axis.function.ref_coord is not None
        ):
            temp = self._lookup.max_exputend(self._storage_name.file_uri)
            if temp:
                chunk.time.axis.function.ref_coord.val = temp
            else:
                self._logger.warning(f'Missing EXPUTEND. Set time to None for {self._storage_name.file_uri}.')
                chunk.time = None
                chunk.time_axis = None


class GHOSTSpatialSpectralTemporal(GHOSTSpectralTemporal):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))
        for extension in range(0, len(self._headers)):
            camera = self._headers[extension].get('CAMERA', extension)
            if camera and camera in ['RED', 'BLUE']:
                # if camera is SLITV, there's no interesting science information in the HDU
                self._accumulate_chunk_position_axes_blueprint(bp, extension)
        self._logger.debug(f'End accumulate_blueprint')


class Gmos(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Observation.instrument.keywords', 'get_provenance_keywords()')
        bp.configure_position_axes((1, 2))

    def get_exposure(self, ext):
        if self.get_obs_type(ext) == 'MASK':
            result = 0.0
        else:
            result = super().get_exposure(ext)
        return result

    def get_time_delta(self, ext):
        if self.get_obs_type(ext) == 'MASK':
            # DB - 05-03-19 - delta hardcoded to 0
            exptime = 0.0
        else:
            exptime = super().get_time_delta(ext)
        return exptime

    def _reset_position(self, observation_type, artifact_type):
        # DB - 04-03-19
        # Another type of GMOS-N/S dataset to archive.
        # Mask images.   json observation_type = “MASK”.
        # These have no WCS info at all, although I guess
        # json ut_date_time could be used as the start date
        # with null exposure time. These would have only
        # instrument, obstype, datatype (spectrum) and
        # product type (AUXILIARY) set.
        return super()._reset_position(observation_type, artifact_type) or observation_type == 'MASK'

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')

        GMOS_RESOLVING_POWER = {
            'B1200': 3744.0,
            'R831': 4396.0,
            'B600': 1688.0,
            'R600': 3744.0,
            'R400': 1918.0,
            'R150': 631.0,
            # current best value available from https://gemini.edu/instrumentation/gmos/components#Gratings
            'B480': 1522.0,  # WF - 16-09-22  ## TODO WF - 23-09-22 - need more accurate min/max values
        }
        # DB 02-12-19
        # Found basic info for the Lya395 filter in a published paper.
        # Central wavelength is .3955 microns and FWHM is 0.00327 microns.
        #
        # The Hartmann ‘filter’ is actually a mask that blocks off half of the
        # light path to the detector but is ‘open’ on the other half.  So same
        # ‘bandpass’ as ‘open’ or ‘empty’.    It’s used when focusing the
        # camera.
        #
        # Blocking filters are only used when doing multi-object
        # spectroscopy.  By reducing the wavelength coverage they are able to
        # add a second row of slits to the mask since the spectra then don’t
        # cover the entire width of the detector in the dispersion direction.

        # 0 = min
        # 1 = max
        # units are microns
        lookup = {
            'GG455': [0.46000, 1.10000],
            'OG515': [0.52000, 1.10000],
            'RG610': [0.61500, 1.10000],
            'RG780': [0.7800, 1.10000],
            'Lya395': [0.393865, 0.397135],
            'open': [0.35000, 1.10000],
            'open2-8': [0.35000, 1.10000],
            'empty': [0.35000, 1.10000],
            'HartmannA': [0.35000, 1.1000],
            'HartmannB': [0.35000, 1.1000],
            'HartmannC': [0.35000, 1.1000],
            'HartmannD': [0.35000, 1.1000],
        }

        # DB - 04-04-19

        # GG455&g filter is a case where both the GG455 order-sorting filter
        # and the g filter are combined (which is likely quite common).  The
        # former is in the ‘lookup’ on line 2343.  The g filter is g_G0301
        # which is in the SVO.  My original ‘svo’ script handled the
        # order-sorting filters as part of that subroutine.  I don’t think
        # that works with the current _update_chunk_energy_gmos.  The o-s
        # filter will change the effective wavelength coverage of the g
        # filter.

        # this means some filter metadata will come from SVO, and some is
        # hard-coded here

        # DB 18-04-19
        #
        # Lots of missing header data for the GMOS observation
        # GN-2005B-Q-60-11-011 so ignore energy WCS

        # DB 02-12-19
        #
        # One problem is 186 files with filter values of empty_0[134].  The
        # ‘empty*’ component should be treated like ‘open’ as far as bandpass is
        # concerned.   And ‘empty_*’ could be stripped out of the displayed filter
        # name. Likely less confusing to remove ‘empty*’ completely.

        reset_energy = False

        filter_md = None
        if 'open' not in filter_name and 'No_Value' not in filter_name and 'empty' not in filter_name:
            filter_md = self._filter_cache.get_filter_metadata(
                self._instrument,
                filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )

        if 'empty' in filter_name:
            # set to 'empty' string here, so can still use lookup logic
            filter_name = re.sub('empty_\\d*', 'empty', filter_name)

        w_max = 10.0
        w_min = 0.0
        for ii in filter_name.split('+'):
            if ii in lookup:
                wl_max = lookup[ii][1]
                wl_min = lookup[ii][0]
                if wl_max < w_max:
                    w_max = wl_max
                if wl_min > w_min:
                    w_min = wl_min

        if filter_md is not None:
            # mingle the SVO lookup values with the hard-coded lookup values
            # from here
            svo_min = (2 * filter_md.central_wl - filter_md.bandpass) / 2
            svo_max = filter_md.bandpass + svo_min
            if svo_max < w_max:
                w_max = svo_max
            if svo_min > w_min:
                w_min = svo_min

        filter_md = svofps.FilterMetadata()
        filter_md.set_bandpass(w_max, w_min)
        filter_md.set_central_wl(w_max, w_min)
        filter_md.set_resolving_power(w_max, w_min)

        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS spectroscopy for {self._storage_name.obs_id}.')
            if math.isclose(filter_md.central_wl, 0.0):
                self._logger.info(
                    f'no spectral wcs, central wavelength is '
                    f'{filter_md.central_wl} for {self._storage_name.obs_id}'
                )
                return
            disperser = self._lookup.disperser(self._storage_name.file_uri)
            # 'unknown' in disperser test obs is GN-2004B-Q-30-15-002
            # 'OLDMIRROR' in disperser test obs is GN-CAL20020329-2-025
            # 'No Value' in disperser test obs is GN-2013B-SV-152-120-001
            if disperser is None or 'unknown' in disperser or 'OLDMIRROR' in disperser or 'No Value' in disperser:
                # DB 23-04-19
                # GMOS observation with OLDMIRROR looks to be a direct image
                # with the pinhole mask in the beam.  But Gemini json info is
                # incorrect since it claims it’s a spectrum.  No energy.
                #
                # GMOS observation GN-CAL20020501-3-000 with unknown
                # disperser.  No energy.  Just the observation number 000
                # flags this one as unusual.
                #
                # DB 13-05-19
                # “No Value” for filter?  Yes, no energy WCS.  Lots of other
                # bogus values in header as well.
                self._logger.warning(f'disperser is {disperser}, no energy.')
                reset_energy = True
            else:
                self.fm = svofps.FilterMetadata()
                self.fm.central_wl = filter_md.central_wl
                self.fm.bandpass = filter_md.bandpass
                if disperser == 'B12000':
                    # DB 16-04-19
                    # B12000 must be a Gemini typo since observation
                    # GN-2006B-Q-39-100-003 has B1200 for the disperser.
                    disperser = 'B1200'
                elif disperser in ['B600-', 'B600+-G5323', '\'B600']:
                    # DB 23-04-19
                    # GMOS observation GS-CAL20030130-1-002 with B600- grating.
                    # The ‘-’ must be a typo.  The observation preceding this
                    # one, GS-CAL20030130-1-001, has the B600 grating.
                    # DB 13-05-19
                    # GS-2013A-Q-91-194-003 S20130517S0092
                    # Yes, disperser value should be B600.  Typo using +-
                    # instead of the usual underscore when someone entered the
                    # info in a config file perhaps?
                    disperser = 'B600'
                if disperser in GMOS_RESOLVING_POWER:
                    self.fm.resolving_power = GMOS_RESOLVING_POWER[disperser]
                else:
                    self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                    raise mc.CadcException(
                        f'{self._instrument}: mystery disperser {disperser} for {self._storage_name.obs_id}'
                    )
        elif data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging for {self._storage_name.obs_id}.')
            self.fm = filter_md
        else:
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        if reset_energy:
            cc.reset_energy(chunk)
        else:
            self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')


class Gnirs(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))

    def get_data_product_type(self, ext):
        # DB - 03-04-19
        # if disperser for GNIRS = MIRROR
        # then it’s an image, otherwise a spectrum.
        result = DataProductType.SPECTRUM
        if self._lookup.disperser(self._storage_name.file_uri) == 'MIRROR':
            result = DataProductType.IMAGE
        return result

    def _reset_energy(self, observation_type, filter_name):
        # DB 16-04-19
        # Energy WCS should be ignored for ‘Moving’ since we don’t know
        # what might have been in the light path when the exposure was
        # actually being taken.
        #
        # 'Dark' test obs is GN-2013B-Q-93-147-036
        # DB 13-05-19
        # GNIRS “Dark” trumps the “L” filter, so no energy.
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and ('Moving' in filter_name or 'Dark' in filter_name or 'DARK' in filter_name)
        )

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin update_energy')

        # DB - 07-02-19
        # Spectroscopy:
        #
        # a) long-slit
        #
        # (Note: spatial WCS info is in the primary header for GNIRS apparently)
        #
        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy is the
        # relevant web page.  Grating, Order-Band, Blocking Filter Range and the
        # two Resolving Power columns (one for each camera) are the important
        # columns.
        #
        # Need to know:
        #
        # grating: from json ‘disperser’ (need 10, 32 or 111 numbers only)
        # Note: if disperser contains string “XD” then the spectrum is cross
        # dispersed.  See b) below.
        #
        # filter: from json ‘filter_name’
        # camera: Short or Long substring from json ‘camera’ value (i.e. “Blue”
        # or “Red” aren’t important)
        # central wavelength: from json ‘central_wavelength’
        # NAXIS2 extension header value
        #
        # crval = central_wavelength
        #
        # use ‘Blocking Filter Range’ for the appropriate filter
        # (e.g. X, J, H...) to determine min/max wavelengths as the rough
        # estimate of the wavelength coverage.

        # b) cross-dispersed (when ‘disperser’ contains ‘?XD’ string)
        #
        # Wavelength ranges given here:
        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/crossdispersed-xd-spectroscopy/xd-prisms
        #
        # You need the ‘Short’ or ‘Long’ from the json ‘camera’ value and the
        # ‘SXD’ or ‘LXD’ from the ‘disperser’ value to look up the coverage.
        #
        # Imaging:
        #
        # Not done very often.  Filter info is here:
        # https://www.gemini.edu/sciops/instruments/gnirs/imaging.  For filter
        # ID’s I think json ‘filter_name’ value of J_(order_5) corresponds to
        # “J (order blocking)” in this table.  ‘filter_name’ value of J
        # corresponds to “J (Mauna Kea)“.

        # DB - 01-03-19
        # Resolution is in the last two columns of the table here, a different
        # value for each grating/filter combination as well as camera value:
        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy.
        # Add that resolution for each grating/filter (and camera) combination
        # sort of as \for the bandpasses.  Already know the grating value.
        # Then, for both long-slit and cross-dispersed spectroscopy need to look
        # at json ‘camera’ (or header CAMERA) value (does it contain ‘short’ or
        # ‘long’) and the json focal_plane_mask value (or SLIT header value;
        # you need only the numeric value at the start of the string).  The
        # camera value tells you which column to look under for the resolution.
        # The best estimate of the resolution, R, is then given by:
        #
        # R = 0.3 x tabulated value / slit width
        #       (for ‘short’ camera observations)
        # R = 0.1 x tabulated value / slit width
        #       (for ‘long’ camera observations)
        #
        # i.e. the resolution gets lower when a wider slit is in place
        # (the tabulated values are for 0.3" and 0.1" slits for the short/long
        # cameras respectively, hence the numbers in the start of each
        # equation). Basically the camera is producing an image of the slit at
        # each wavelength but if the slit is wider the dispersed image is also
        # wider and so different colours are blended together more so the
        # resolution gets worse as you open up the slit (but you get more light
        # through a wider slit so it’s a compromise between throughput and
        # spectral resolution).

        # cross-dispersion information:
        # Change xd_mode to include grating ID (e.g. ‘32’ or ‘10’ or ‘110’)
        # and two more configurations:
        # xd_mode = {‘SB+SXD+32’: [0.9, 2.5, 1800.0],
        # ‘LB+LXD+10’: [0.9, 2.5, 1800.0],
        # ‘LB+SXD+10’: [1.2, 2.5, 1800.0],
        # ‘SB+SXD+110’:[0.9, 2.5, 5400.0],
        # ‘LB+SXD+110’:[0.9, 2.5, 5400.0],
        # ‘SB+LXD+32’:{0.9, 2.5, 5400.0]}

        # DB - 05-03-19
        # Disregard central wavelength provided by Gemini.  e.g. if a long-slit
        # observation with the K filter in the beam has a central_wavelength
        # setting of (for example) 2.3 microns then using the bounding
        # wavelengths in the arrays to define the wavelength coverage will
        # result in the upper wavelength being outside the filter bandpass.
        # Instead use the average of the lower/upper wavelength ranges
        # for each configuration in long_slit_mode and xd_mode as the central
        # wavelength.  e.g. for K long-slit observations the central wavelength
        # would be (1.91+2.49)/2.0 or 2.2 micro^ns.  The way it is now with a
        # central wavelength of 2.3 microns we calculate wavelength limits of
        # 2.01 to 2.59 microns but wavelengths beyond 2.49 microns don’t make it
        # past the filter.

        # DB - 03-04-19

        # GNIRS observation GN-2017A-Q-44-25-031 has a json ‘mode’ value of
        # ‘imaging’ despite having a disperser value of ‘32_mm&SXD’ indicating
        # that it’s a spectrum (as you can also see from the preview). Means
        # forgetting about relying on json ‘mode’ but looking to see if the
        # json disperser value = ‘MIRROR’.

        # DB 08-09-19 - observations like GN-2010B-SV-142-655-012 with null
        # 'disperser'.  Some kind of rare calibration exposure - there are
        # only on the order of 100 in the archive - or just an error in the
        # header.  The spectra look ‘dispersed’ to me.  Just ignore energy
        # WCS for these.

        # DB 24-04-19
        # ND = neutral density and so any ND* filter can be ignored as it
        # shouldn’t affect transmission band.  Likely observing a bright
        # target and they need to reduce light by a factor of 100X.

        # DB 10-07-19
        # Info is available here:
        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters.
        # This filter is used for acquisition of targets in spectroscopy mode
        # sometimes.  Use the central wavelength of 3.295 microns and width of
        # 1.5% of that or 0.049425 microns.

        # long slit mode information source:
        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy
        # 0 - min wavelength
        # 1 - max wavelength
        # 2 - 'short' camera resolution
        # 3 - 'long' camera resolution
        # 4 - Since November 2012 and for the cross-dispersed mode with
        # the 2 pix wide slit only resolving powers are somewhat lower, as
        # follows: X-1400; J-1400, H-1400; K-1300, for 'short' camera
        # resolution

        # DB 21-08-19
        # GN-2010B-SV-142-300-002:  add a line in the ‘10’ section of
        # gnirs_lookup similar to the ‘L’ entry but with the PAH filter
        # bandpass:  ‘PAH’: [3.270, 3.320, 570, 1800].  These are likely
        # mis-identified acquisition observations rather than science.

        # DB 31-03-21
        # It looks like Gemini changed FITS header content for gratings, etc. at
        # some point in time.  There are two prisms used in cross-dispersed mode
        # and these are identified by SXD and LXD but the leading letter.   If
        # you find an “&XD” entry in the JSON disperser value, then look at the
        # value of the FITS header PRISM = 'SXD_G5509' to find the letter that
        # should be in front of the ‘XD’ in the JSON result? It might not be the
        # first letter in the value: new files have PRISM   = 'SB+SXD_G5536'

        gnirs_lookup = {
            '10': {
                'X': [1.03, 1.17, 570, 2100],
                'J': [1.17, 1.37, 570, 1600],
                'H': [1.47, 1.80, 570, 1700],
                'K': [1.91, 2.49, 570, 1700],
                'L': [2.80, 4.20, 570, 1800],
                'M': [4.40, 6.00, 570, 1200],
                'LB+LXD': [0.9, 2.5, 1800],
                'LB+SXD': [1.2, 2.5, 1800],
                'PAH': [3.270, 3.320, 570, 1800],
            },
            '32': {
                'X': [1.03, 1.17, 1700, 5100, 1400],
                'J': [1.17, 1.37, 1600, 4800, 1400],
                'H': [1.49, 1.80, 1700, 5100, 1400],
                'K': [1.91, 2.49, 1700, 5100, 1300],
                'L': [2.80, 4.20, 1800, 5400, 1800],
                'M': [4.40, 6.00, 1240, 3700, 1240],
                'SB+SXD': [0.9, 2.5, 1800, 1800],
                'SB+LXD': [0.9, 2.5, 5400, 5400],
                'LB+LXD': [0.9, 2.5, 5400, 5400],
                'LB+SXD': [1.2, 2.5, 1800, 1800],
                'PAH': [3.270, 3.320, None, None],
            },
            '111': {
                'X': [1.03, 1.17, 6600, 17800],
                'J': [1.17, 1.37, 7200, 17000],
                'H': [1.49, 1.80, 5900, 17800],
                'K': [1.91, 2.49, 5900, 17800],
                'L': [2.80, 4.20, 6400, 19000],
                'M': [4.40, 6.00, 4300, 12800],
                'SB+SXD': [0.9, 2.5, 5400],
                'LB+SXD': [0.9, 2.5, 5400],
                'LB+LXD': [0.9, 2.5, 17000],
                'PAH': [3.270, 3.320, None, None],
            },
        }

        reset_energy = False
        if 'open' in filter_name:
            # DB 21-08-19
            # GS-2004A-DD-8-10-001:  set the bandpass to between 1 micron and
            # 6 microns, same as for GMOS imaging with ‘open’ filter.  This is
            # just an acquisition observation.
            w_min = 1.0
            w_max = 6.0
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(w_max, w_min)
            self.fm.set_central_wl(w_max, w_min)
            self.fm.set_resolving_power(w_max, w_min)
        else:
            self.fm = svofps.FilterMetadata('GNIRS')
            if data_product_type == DataProductType.SPECTRUM:
                self._logger.info(f'SpectralWCS Spectroscopy mode for {self._storage_name.obs_id}.')
                disperser = self._lookup.disperser(self._storage_name.file_uri)
                if disperser is None:
                    self._logger.info(f'No disperser. No energy for {self._storage_name.obs_id}')
                    reset_energy = True
                else:
                    grating = disperser.split('_')[0]
                    # 'UNKNOWN' in grating test obs GS-CAL20040924-6-006
                    # 'ENG -' in grating test obs GN-CAL20130813-22-010
                    if 'UNKNOWN' in grating or 'Moving' in grating or 'ENG' in grating:
                        # DB 23-04-19 - GNIRS grating UNKNOWN:  no energy
                        # DB 18-04-19
                        # If grating is moving then the observation is almost
                        # certainly garbage so ignore energy WCS.
                        # DB 13-05-19
                        # No energy for the invalid “ENG - 170000" grating
                        # entry.
                        self._logger.warning(f'grating is {grating}. No energy for {self._storage_name.obs_id}')
                        reset_energy = True
                    else:
                        if grating not in gnirs_lookup:
                            self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                            raise mc.CadcException(
                                f'{self._instrument} Mystery grating {grating} for {self._storage_name.obs_id}'
                            )

                        camera = self._lookup.camera(self._storage_name.file_uri)
                        if camera == 'Moving':
                            # DB 21-08-19
                            # No energy for the ‘moving’ camera
                            self._logger.info(f'Camera is moving. No energy for {self._storage_name.obs_id}')
                            reset_energy = True
                        else:
                            focal_plane_mask = self._lookup.focal_plane_mask(self._storage_name.file_uri)
                            if focal_plane_mask is None or 'arcsec' not in focal_plane_mask:
                                # DB 24-04-19
                                # Assume slit width of 1.0 for GNIRS
                                # observation without a focal plane mask
                                # value.
                                slit_width = 1.0
                            else:
                                slit_width = mc.to_float(focal_plane_mask.split('arcsec')[0])

                            if 'XD' in disperser:
                                self._logger.debug(f'cross dispersed mode for {self._storage_name.obs_id}.')
                                # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/crossdispersed-xd-spectroscopy/xd-prisms
                                #
                                # 0 = lower
                                # 1 = upper
                                # 2 = spectral resolution with 2-pix wide slit
                                # DB - 04-03-19 - Change the last number in
                                # each row to 1800.0 since the resolving power
                                # is the same for all 3 cases
                                #
                                # Add line to find grating ID (from long-slit
                                # code):
                                # grating = disperser.split(‘_’)[0]
                                #
                                # Then change ‘lookup’ to include grating.
                                #
                                # I can’t find any other combinations (e.g.
                                # ‘LB+LXD+32’) but no guarantee that I won’t
                                # have to add another line or two if we see
                                # failures. Wavelength coverage isn’t correct
                                # for the R=5400 entries because only about
                                # 1/3rd of the full band pass is covered but
                                # in bits and pieces that we can’t identify in
                                # raw image.

                                lookup = None
                                if '&XD' in disperser:
                                    prism = self._headers[0].get('PRISM')
                                    cd_mode = prism[prism.index('XD') - 1]
                                    coverage = f'{cd_mode}XD'
                                else:
                                    coverage = disperser[-3:]

                                if camera.startswith('Short'):
                                    lookup = f'SB+{coverage}'
                                elif camera.startswith('Long'):
                                    lookup = f'LB+{coverage}'

                                if camera.startswith('Long'):
                                    slit_table_value = 0.1
                                    lookup_index = 2
                                elif camera.startswith('Short'):
                                    slit_table_value = 0.3
                                    lookup_index = 2
                                else:
                                    self._observable.rejected.record(
                                        mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name
                                    )
                                    raise mc.CadcException(
                                        f'{self._instrument} Mystery camera '
                                        f'definition {camera} for '
                                        f'{self._storage_name.obs_id}'
                                    )
                            else:
                                self._logger.debug(f'long slit mode for {self._storage_name.obs_id}.')
                                if filter_name == 'PAH':
                                    lookup = filter_name
                                else:
                                    bandpass = filter_name[0]
                                    lookup = bandpass
                                    if camera.startswith('Long'):
                                        slit_table_value = 0.1
                                        lookup_index = 3
                                    elif camera.startswith('Short'):
                                        date_time = ac.get_datetime_mjd(
                                            self._lookup.ut_datetime(self._storage_name.file_uri)
                                        )
                                        if date_time > ac.get_datetime_mjd('2012-11-01T00:00:00'):
                                            slit_table_value = 0.3
                                            if grating == '32':
                                                lookup_index = 4
                                            else:
                                                lookup_index = 3
                                        else:
                                            slit_table_value = 0.3
                                            lookup_index = 2
                                    else:
                                        self._observable.rejected.record(
                                            mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name
                                        )
                                        raise mc.CadcException(
                                            f'{self._instrument}: Mystery '
                                            f'camera definition {camera} for '
                                            f'{self._storage_name.obs_id}'
                                        )

                            lookup = lookup.upper()
                            if lookup not in gnirs_lookup[grating]:
                                self._observable.rejected.record(
                                    mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name
                                )
                                raise mc.CadcException(
                                    f'{self._instrument}: Mystery lookup '
                                    f'{lookup} for grating {grating}, obs '
                                    f'{self._storage_name.obs_id}'
                                )
                            bounds = gnirs_lookup[grating][lookup]
                            self.fm.set_bandpass(bounds[1], bounds[0])
                            if lookup == 'PAH':
                                self.fm.resolving_power = None
                            else:
                                self.fm.resolving_power = slit_table_value * bounds[lookup_index] / slit_width
                            self.fm.set_central_wl(bounds[1], bounds[0])
            elif data_product_type == DataProductType.IMAGE:
                self._logger.debug(f'SpectralWCS imaging mode for {self._storage_name.obs_id}.')
                # https://www.gemini.edu/sciops/instruments/gnirs/imaging
                # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters

                # DB 23-04-19
                # GNIRS filter “L”: add lower/upper bandpass info for the
                # L and M filter on this page,
                # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters
                # to ‘imaging’ dictionary.  (In ‘imaging’ I think the K order
                # blocking filter lower bandpass should be 1.91).
                # Despite Gemini’s footnote stating that L and M filters are
                # NOT used for acquisition images.

                # 0 - min wavelength
                # 1 - max wavelength
                imaging = {
                    'Y': [0.97, 1.07],
                    'J': [1.17, 1.33],
                    'J order blocking)': [1.17, 1.37],
                    'H': [1.49, 1.80],
                    'K': [2.03, 2.37],
                    'K order blocking': [1.91, 2.49],
                    'H2': [2.105, 2.137],
                    'PAH': [3.27, 3.32],
                    'X': [1.03, 1.17],
                    'XD': [0.9, 2.56],
                    'M': [4.4, 6.0],
                    'L': [2.8, 4.2],
                }

                # DB 30-04-19
                # GN-CAL20180607-27-001:  ignore ‘PupilViewer’ if possible.
                # DB 21-08-19
                # PV+K: ‘PV’ = pupil viewer.  Shouldn’t affect transmission
                # values (I think) so ignore when looking up lower/upper
                # bandpasses.
                filter_name = filter_name.replace('PupilViewer', '').replace('PV', '')
                filter_name = filter_name.strip('+')

                if len(filter_name) == 1 or filter_name == 'H2' or filter_name == 'PAH' or filter_name == 'XD':
                    bandpass = filter_name
                else:
                    bandpass = filter_name[0]

                bandpass = bandpass.upper()
                if bandpass in imaging:
                    bounds = imaging[bandpass]
                else:
                    raise mc.CadcException(
                        f'{self._instrument}: Unexpected filter_name '
                        f'{filter_name} for {self._storage_name.obs_id}'
                    )

                self.fm.set_central_wl(bounds[1], bounds[0])
                self.fm.set_bandpass(bounds[1], bounds[0])
            else:
                raise mc.CadcException(
                    f'{self._instrument}: Unexpected DataProductType '
                    f'{data_product_type} for {self._storage_name.obs_id}'
                )

        if reset_energy:
            cc.reset_energy(chunk)
        else:
            self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')


class Gpi(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        equinox = self._headers[0].get('TRKEQUIN')
        if equinox and 1800.0 <= equinox <= 2500.0:
            # DB 07-06-21
            # No spatial WCS if equinox not properly set
            #
            # DB - 22-02-19 - GPI
            # WCS info is garbage in header, so for spatial WCS (for both image and spectrum)
            # override with the following functions
            bp.configure_position_axes((1, 2))
            for extension in range(1, len(self._headers)):
                bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN', extension)
                bp.set('Chunk.position.axis.axis1.cunit', 'deg', extension)
                bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN', extension)
                bp.set('Chunk.position.axis.axis2.cunit', 'deg', extension)
                bp.set('Chunk.position.axis.function.cd11', 'get_cd11()', extension)
                # cd1_2 = cd2_1 = 0.0 since we don’t know rotation value
                bp.set('Chunk.position.axis.function.cd12', 0.0, extension)
                bp.set('Chunk.position.axis.function.cd21', 0.0, extension)
                bp.set('Chunk.position.axis.function.cd22', 'get_cd22()', extension)
                bp.set('Chunk.position.axis.function.refCoord.coord1.pix', 'get_crpix1()', extension)
                # crval1 = json ra value
                # crval2 = json dec value
                bp.set('Chunk.position.axis.function.refCoord.coord1.val', 'get_ra()', extension)
                bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 'get_crpix2()', extension)
                bp.set('Chunk.position.axis.function.refCoord.coord2.val', 'get_dec()', extension)
                bp.set('Chunk.position.equinox', '_get_trkequin()', extension)
                bp.add_attribute('Chunk.position.coordsys', 'RADESYS', extension)

    def _get_trkequin(self, ext):
        return self._headers[0].get('TRKEQUIN')

    def get_cd11(self, ext, keyword='NAXIS1'):
        # FOV is 2.8" x 2.8" on each side or 2.8/3600.0 degrees
        # cd1_1 = 2.8/(3600.0 * naxis1)
        # cd2_2 = 2.8/(3600.0 * naxis2)
        return RADIUS_LOOKUP[self._instrument] / self._headers[ext].get(keyword)

    def get_cd22(self, ext):
        return self.get_cd11(ext, 'NAXIS2')

    def get_crpix1(self, ext, keyword='NAXIS1'):
        # naxis1 = extension header NAXIS1 value
        # naxis2 = extension header NAXIS2 value
        # crpix1 = naxis1/2.0
        # crpix2 = naxis2/2.0
        naxis1 = self._headers[ext].get(keyword)
        if naxis1:
            result = naxis1 / 2.0
        return result

    def get_crpix2(self, ext):
        return self.get_crpix1(ext, 'NAXIS2')

    def get_data_product_type(self, ext):
        mode = self._lookup.mode(self._storage_name.file_uri)
        if mode is None:
            raise mc.CadcException(
                f'{self._instrument}: No mode information found for {self._storage_name.file_name}'
            )
        # DB - 22-02-19 FOR GPI only:  To determine if the data type
        # is an ‘image’ or ‘spectrum’:
        #     json ‘mode’ = IFP then ‘image’
        #     json ‘mode’ = IFS then ‘spectrum’
        if mode in ['IFP', 'imaging']:
            result = DataProductType.IMAGE
        elif mode == 'IFS':
            result = DataProductType.SPECTRUM
        else:
            self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
            raise mc.CadcException(f'{self._instrument} Mystery mode {mode} for {self._storage_name.file_name}')
        return result

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')

        # DB - 22-02-19
        # For imaging energy WCS, use standard imaging algorithm. i.e central
        # wavelength, bandpass from SVO filter, and
        # resolution = central_wavelength/bandpass
        #
        # For energy WCS:
        #
        # naxis = extension header NAXIS2 value
        #
        # Need to hard-code the resolution from values in top table on this
        # page:
        # https://www.gemini.edu/sciops/instruments/gpi/observing-modes-metamodes.
        # Use mean of column 4 range.  e.g. for H filter use 46.5
        # units are microns
        gpi_lookup = {
            'Y': (36 + 34) / 2.0,
            'J': (39 + 35) / 2.0,
            'H': (49 + 44) / 2.0,
            'K1': (70 + 62) / 2.0,
            'K2': (83 + 75) / 2.0,
        }

        filter_md = self._filter_cache.get_filter_metadata(
            self._instrument,
            filter_name,
            self._lookup.telescope(self._storage_name.file_uri),
        )
        if data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging mode for {self._storage_name.obs_id}.')
            self.fm = filter_md
        elif data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS Spectroscopy mode for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = filter_md.central_wl
            self.fm.bandpass = filter_md.bandpass
            if filter_name in gpi_lookup:
                self.fm.resolving_power = gpi_lookup[filter_name]
            else:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(
                    f'{self._instrument}:  Mystery filter name {filter_name} '
                    f'for resolving power {self._storage_name.obs_id}'
                )
        else:
            raise mc.CadcException(
                f'{self._instrument}:  Do not understand DataProductType '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )

        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_position(self, part, chunk, extension):
        pass


class GracesSpectralTemporal(GeminiMapping):
    """Fibre-fed spectrograph"""

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.set('Plane.provenance.lastExecuted', 'get_provenance_last_executed()')
        bp.set('Plane.provenance.producer', 'get_provenance_producer()')
        bp.set('Plane.provenance.reference', 'get_provenance_reference()')
        bp.set('Plane.provenance.version', 'get_provenance_version()')

    def _get_filter_name(self, ext):
        return None

    def _reset_position(self, observation_type, artifact_type):
        pass

    def _update_energy(self, chunk, data_product_type, filter_name):
        """graces-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin _update_energy')

        # DB - 22-02-19  Axis 2 is the spectral axis.   Resolution is
        # 67,500 divided by second number in json ‘detector_binning’
        # value of ‘N x N’.  The central wavelength in json
        # is always set to 0.7 microns and range is from about 0.4 to
        # 1.0 microns.  So use 0.7 as the crval, crpix = naxis2/2.0 and
        # delta = (1.0 - 0.4)/naxis2.  Just a kludge to get spectral
        # range correct for an echelle spectrograph.

        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS {data_product_type} for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = self._lookup.central_wavelength(self._storage_name.file_uri)
            self.fm.set_bandpass(1.0, 0.4)
            ccd_sum = self._lookup.detector_binning(self._storage_name.file_uri)
            self.fm.resolving_power = 67500.0
            if ccd_sum is not None:
                temp = float(ccd_sum.split('x')[1])
                self.fm.resolving_power = 67500.0 / temp
            self._build_chunk_energy(chunk, filter_name)
        else:
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._logger.debug(f'End _update_energy')

    def _update_position(self, part, chunk, ext):
        pass

    def _update_time(self, chunk):
        cc.undo_astropy_cdfix_call(chunk, self.get_time_delta(0))


class GracesSpatialSpectralTemporal(GracesSpectralTemporal):
    """Fibre-fed spectrograph"""

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))
        self._accumulate_chunk_position_axes_blueprint(bp, 0)


class Gsaoi(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))

    def _reset_energy(self, observation_type, filter_name):
        # DB 24-04-19
        # ‘Unknown+Blocked2’ filter, no spectral WCS.
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and ('Unknown' in filter_name or 'Blocked' in filter_name)
        )

    def _update_energy(self, chunk, data_product_type, filter_name):
        """General chunk-level Energy WCS construction.
        Same as Fox.
        """
        self._logger.debug(f'Begin _update_energy')
        if data_product_type is DataProductType.IMAGE:
            self._logger.debug(f'Spectral WCS {data_product_type} mode for {self._storage_name.obs_id}.')
            self.fm = self._filter_cache.get_filter_metadata(
                self._instrument,
                filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )
            if self.fm is None:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(f'{self._instrument}: mystery filter {filter_name}')
            self._build_chunk_energy(chunk, filter_name)
        else:
            raise mc.CadcException(
                f'{self._instrument} no Spectral WCS support when '
                f'DataProductType {data_product_type} for '
                f'{self._storage_name.obs_id}'
            )
        self._logger.debug(f'End _update_energy')


class Hokupaa(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')

    def get_art_product_type(self, ext):
        result = ProductType.SCIENCE
        image_type = self._headers[ext].get('IMAGETYP')
        # took the list from the AS GEMINI Obs. Type values
        if image_type in ['DARK', 'BIAS', 'CAL', 'ARC', 'FLAT']:
            result = ProductType.CALIBRATION
        return result

    def get_cd11(self, ext):
        pix_scale = self._headers[ext].get('PIXSCALE')
        if pix_scale is None:
            result = None
        else:
            result = mc.to_float(pix_scale) / 3600.0
        return result

    def get_cd22(self, ext):
        return self.get_cd11(ext)

    def get_crpix1(self, ext, keyword='CRPIX1'):
        naxis1 = self._headers[ext].get(keyword)
        if naxis1 is None:
            result = None
        else:
            result = naxis1 / 2.0
        return result

    def get_crpix2(self, ext):
        return self.get_crpix1(ext, 'CRPIX2')

    def get_dec(self, ext):
        ra_ignore, dec = self.get_sky_coord(ext, 'RA', 'DEC')
        return dec

    def get_obs_type(self, ext):
        # DB - 01-18-19 Use IMAGETYP as the keyword
        return self._headers[ext].get('IMAGETYP')

    def get_ra(self, ext):
        ra, dec_ignore = self.get_sky_coord(ext, 'RA', 'DEC')
        return ra

    def _reset_energy(self, observation_type, filter_name):
        # DB 14-08-19
        # no energy since ‘home’ might either be an open position or a blocked
        # position.
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and ('LowFlx' in filter_name or 'home' in filter_name)
        )

    def _reset_position(self, observation_type, artifact_type):
        result = super()._reset_position(observation_type, artifact_type)
        ra = self.get_ra(0)
        if ra is None:
            result = True
        return result

    def _update_energy(self, chunk, data_product_type, filter_name):
        """hokupaa-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin _update_energy')

        # DB - 18-01-19 - Note: it is always imaging for Hokupa'a + QUIRC so
        # that could be hardcoded.

        # DB - 20-02-19 - This page has a table of QUIRC (the camera part of
        # Hokupa’a + QUIRC) filter names and bandpasses (in microns).
        # Use this as a lookup for central wavelengths and bandpass
        # http://www.gemini.edu/sciops/instruments/uhaos/uhaosQuirc.html

        # J+CO is "Dark Position", units are microns
        # 0 - central wavelength
        # 1 - bandpass
        hokupaa_lookup = {
            'J': [1.25, 0.171],
            'H': [1.65, 0.296],
            'K': [2.2, 0.336],
            'K\'': [2.12, 0.41],
            'Kp': [2.12, 0.41],
            'H+K notch': [1.8, 0.7],
            'methane low': [1.56, 0.0120],
            'methane high': [1.71, 0.0120],
            'FeII': [1.65, 0.0170],
            'HeI': [2.06, 0.0030],
            '1-0 S(1) H2': [2.12, 0.0023],
            'H Br(gamma)': [2.166, 0.0150],
            'K-continuum': [2.26, 0.0060],
            'CO': [2.29, 0.0020],
        }

        # DB 27-05-19
        # The bandpasses for the 7 bottom entries have to be corrected
        # because of the misleading column heading.  The bandpasses for these
        # 7 are in Angstroms, not microns.  So divide the values by 10,000 to
        # convert to microns. 170 Å = 0.017 microns for the FeII filter. H2/23
        # = 1-0 S(1) H2 - bandpass should be 0.0023
        # HKnotch = H+K notch
        # 1.56/120 = methane low  - bandpass should be 0.0120
        # 1.71/120 = methane high - ditto
        # LowFlx is likely the J+CO combo in line 14 of the filter table in
        # the web page:  effectively a shutter so no energy in this case. And I
        # think FeI/17 is likely meant to be the FeII filter (just going by the
        # “17” in the bandpass.  I’m not aware of an IR Fe I filter…

        # DB 21-08-19
        # GN-2001A-C-23-38-613:  filter BrG/20 -> BrG = not sure what the
        # ‘/20’ denotes since it is not the correct bandpass if that’s what
        # it’s for, but likely Br(gamma) filter, #11 on the list on the web
        # page
        #
        # GN-2001B-C-3-75-110:   filter HeI/30 -> He I], #9 on the list.  In
        # this case ‘/30’ might be the bandpass  RA/Dec values in the header
        # look reasonable but archive suggests ‘invalid’ for some reason so no
        # spatial WCS I guess.  DEC value might not agree with the telescope
        # elevation value…
        #
        # GN-2001B-C-16-19-039:  2.26/60 - it’s K-continuum.
        #
        # GN-2001A-DD-2-1-537:  Much of the header info seems to be incorrect
        # so no WCS.  But what’s odd is that searching Gemini’s archive for
        # this data label returns 3 files with different file names:  23, 24
        # and 25 February.  EQUINOX only seems to be incorrect in
        # 01FEB23_537.fits.

        filter_name_repair = {
            'H2/23': '1-0 S(1) H2',
            'HKnotch': 'H+K notch',
            '1.56/120': 'methane low',
            '1.71/120': 'methane high',
            'FeI/17': 'FeII',
            '2.26/60': 'K-continuum',
            'BrG/20': 'H Br(gamma)',
            'HeI/30': 'HeI',
        }

        if filter_name in filter_name_repair:
            filter_name = filter_name_repair[filter_name]
        if filter_name not in hokupaa_lookup:
            self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
            raise mc.CadcException(
                f'{self._instrument}: Mystery filter {filter_name} for {self._storage_name.obs_id}'
            )
        if data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging mode for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = hokupaa_lookup[filter_name][0]
            self.fm.bandpass = hokupaa_lookup[filter_name][1]
        else:
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug(f'End _update_energy')

    def _update_position(self, part, chunk, ext):
        header = self._headers[ext]
        equinox = mc.to_float(header.get('EQUINOX'))
        if not 1800.0 <= equinox <= 2500.0:
            self._logger.warning(
                f'EQUINOX value is wrong ({equinox}), no spatial WCS for {self._storage_name.obs_id}'
            )
            return

        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(ext)
        header['CRVAL2'] = self.get_dec(ext)
        header['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        header['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = self.get_crpix1(ext)
        header['CRPIX2'] = self.get_crpix2(ext)
        header['CD1_1'] = self.get_cd11(ext)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(ext)

        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, ext)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        self._logger.debug('End _update_position')


class Hrwfs(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')
        bp.configure_position_axes((1, 2))

    def _get_filter_name(self, ext):
        filter_name = super()._get_filter_name(ext)
        if filter_name is None or len(filter_name.strip()) == 0:
            filter_name = self._search_through_keys(0, ['FILTER1', 'FILTER2'])
        self._logger.info(f'Filter names are "{filter_name}" in {self._storage_name.obs_id}')
        return filter_name

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        # DB 27-05-19
        # e.g. GS-CAL20020322-7-0003 2002mar22_0055, filter1=‘neutral’,
        # filter2=‘open’, so treat similarly to “open” for GMOS since it’s a
        # similar CCD:  maybe 0.35 - 1.0 microns.
        if data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'Spectral WCS {data_product_type} mode for {self._storage_name.obs_id}.')
            if 'open' in filter_name or 'neutral' in filter_name or 'undefined' in filter_name:
                w_min = 0.35
                w_max = 1.0
                self.fm = svofps.FilterMetadata()
                self.fm.set_bandpass(w_max, w_min)
                self.fm.set_central_wl(w_max, w_min)
                self.fm.set_resolving_power(w_max, w_min)
                # DB 27-05-19
                # bandpassName likely best to set to NULL
                filter_name = ''
            else:
                self.fm = self._filter_cache.get_filter_metadata(
                    self._instrument,
                    filter_name,
                    self._lookup.telescope(self._storage_name.file_uri),
                )
                temp = []
                for ii in filter_name.split('+'):
                    temp.append(ii.split('_')[0])
                filter_name = '+'.join(i for i in temp)
            self._build_chunk_energy(chunk, filter_name)
        else:
            raise mc.CadcException(
                f'{self._instrument} no Spectral WCS support when '
                f'DataProductType {data_product_type} for '
                f'{self._storage_name.obs_id}'
            )
        self._logger.debug('End _update_energy')


class Igrins(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        telescope = self._lookup.telescope(self._storage_name.file_uri)
        if telescope is not None and 'North' in telescope:
            x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
        else:
            x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
        bp.set('Observation.telescope.geoLocationX', x)
        bp.set('Observation.telescope.geoLocationY', y)
        bp.set('Observation.telescope.geoLocationZ', z)
        # avoid the self._lookup.data_label, because that data label is
        # repaired with DDMM to be more unique
        data_label = self._headers[0].get('DATALAB')
        bp.set(
            'Plane.provenance.reference',
            f'http://archive.gemini.edu/searchform/{data_label}',
        )
        bp.configure_position_axes((1, 2))

    def get_target_type(self, ext):
        return TargetType.OBJECT

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        # 0 minimum wavelength in microns
        # 1 max wavelength in microns
        lookup = {'H': [1.47, 1.83], 'K': [1.93, 2.5]}
        if filter_name in lookup.keys():
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(lookup.get(filter_name)[1], lookup.get(filter_name)[0])
            self.fm.set_central_wl(lookup.get(filter_name)[1], lookup.get(filter_name)[0])
            self.fm.resolving_power = 45000.0
        else:
            self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
            raise mc.CadcException(
                f'{self._instrument}: Mystery filter {filter_name} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_position(self, part, chunk, ext):
        self._logger.debug('Begin _update_position')
        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(0)
        header['CRVAL2'] = self.get_dec(0)
        header['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        header['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = 1.0
        header['CRPIX2'] = 1.0
        header['CD1_1'] = RADIUS_LOOKUP[self._instrument]
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = RADIUS_LOOKUP[self._instrument]

        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, ext)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        self._logger.debug('End _update_position')


class Igrins2(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        telescope = self._lookup.telescope(self._storage_name.file_uri)
        if telescope is not None and 'North' in telescope:
            x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
        else:
            x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
        bp.set('Observation.telescope.geoLocationX', x)
        bp.set('Observation.telescope.geoLocationY', y)
        bp.set('Observation.telescope.geoLocationZ', z)
        # avoid the self._lookup.data_label, because that data label is
        # repaired with DDMM to be more unique
        data_label = self._headers[0].get('DATALAB')
        bp.set('Plane.provenance.reference', f'http://archive.gemini.edu/searchform/{data_label}')
        bp.configure_position_axes((1, 2))

    def get_target_type(self, ext):
        return TargetType.OBJECT

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        # https://gemini.edu/instrumentation/igrins-2/capability
        # 0 minimum wavelength in microns
        # 1 max wavelength in microns
        lookup = {'H': [1.49, 1.80], 'K': [1.96, 2.46]}
        if filter_name in lookup.keys():
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(lookup.get(filter_name)[1], lookup.get(filter_name)[0])
            self.fm.set_central_wl(lookup.get(filter_name)[1], lookup.get(filter_name)[0])
            self.fm.resolving_power = 45000.0
        else:
            self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
            raise mc.CadcException(
                f'{self._instrument}: Mystery filter {filter_name} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_artifact(self, artifact):
        self._logger.debug(f'Begin _update_artifact for {artifact.uri}')
        for part in artifact.parts.values():
            for chunk in part.chunks:
                chunk.naxis = None
                if chunk.time:
                    chunk.time_axis = None
                if chunk.position:
                    chunk.position_axis_1 = None
                    chunk.position_axis_2 = None
        self._logger.debug('End _update_artifact')

    def _update_position(self, part, chunk, ext):
        self._logger.debug('Begin _update_position')
        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(0)
        header['CRVAL2'] = self.get_dec(0)
        # https://gemini.edu/instrumentation/igrins-2
        header['CDELT1'] = 0.33  # arcseconds
        header['CDELT2'] = 5.0  # arcseconds
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = 1.0
        header['CRPIX2'] = 1.0
        header['CD1_1'] = 0.33  # arcseconds
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = 5.0  # arcseconds

        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, ext)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = None
        chunk.position_axis_2 = None
        self._logger.debug('End _update_position')


class MAROONXTemporal(GeminiMapping):
    """
    A high-resolution (R~80,000) optical (500-920nm), bench-mounted, fiber-fed echelle spectrograph designed to
    deliver 1 m/s radial velocity precision for M dwarfs down to and beyond V = 16.
    https://gemini.edu/instrumentation/maroon-x
    https://gemini.edu/instrumentation/maroon-x/capability
    https://astro.uchicago.edu/~jbean/spectrograph.html
    https://arxiv.org/pdf/1606.07140
    https://astro.uchicago.edu/~jbean/documents/SPIE2018_MAROON_X.pdf

    WF 07-01-25
    It is a spectrograph like any other instrument, and so its observations are fully described by a wavelength
    range, and resolution. It has two arms (blue and red; just like ghost), so that each integration comes along
    with two separate exposures. The resolution is typically R~85k. The data from both arms are contained in each
    MEF image, so the spectral range can be set to 500-920nm.

    Both the spectral range and the resolution are fixed, and will always be the same. Similarly, there are no user
    selectable filters. So in that sense, the data will always be the same, with just times and locations changing
    from image to image.
    """

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        self._accumulate_obs_plane_artifact_blueprint(bp)

        bp.configure_time_axis(4)
        bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.error.syser', '1e-07')
        bp.set('Chunk.time.axis.error.rnder', '1e-07')
        bp.set('Chunk.time.axis.function.naxis', '1')
        bp.set('Chunk.time.axis.function.refCoord.val', 'get_time_function_val()')
        bp.add_attribute('Chunk.time.timesys', 'TIMESYS')
        bp.set_default('Chunk.time.timesys', 'UTC')
        for extension in range(0, len(self._headers)):
            camera_arm = self._headers[extension].get('ARM')
            if camera_arm and camera_arm in ['RED', 'BLUE']:
                bp.set('Chunk.time.resolution', 'get_exposure()', extension)
                bp.set('Chunk.time.exposure', 'get_exposure()', extension)
                bp.set('Chunk.time.axis.function.delta', 'get_time_delta()', extension)
        self._logger.debug(f'End accumulate_blueprint')

    def update(self):
        try:
            cc.TelescopeMapping2.update(self)
            if self._headers is None:
                # proprietary header metadata at archive.gemini.edu, so do nothing
                return self._observation

            # processed files
            if cc.is_composite(self._headers) and not isinstance(self._observation, DerivedObservation):
                self._observation = self._update_composite(self._observation)

            if isinstance(self._observation, DerivedObservation):
                cc.update_observation_members(self._observation)

            if self._observation.proposal is not None and self._observation.proposal.pi_name is None:
                program = self._storage_name.md_context.pi_metadata.get(self._observation.proposal.id)
                if program:
                    self._observation.proposal.pi_name = program.get('pi_name')
                    self._observation.proposal.title = program.get('title')

            GeminiMapping.value_repair.repair(self._observation)
        except Exception as e:
            self._logger.error(f'Error {e} for {self._observation.observation_id}')
            tb = traceback.format_exc()
            self._logger.debug(tb)
            raise mc.CadcException(e)
        self._logger.debug('Done update.')
        return self._observation

    def get_exposure(self, ext):
        return _to_float(self._headers[ext].get('EXPTIME'))

    def get_time_delta(self, ext):
        result = None
        exptime = self.get_exposure(ext)
        if exptime:
            result = mc.convert_to_days(exptime)
        return result

    def _update_plane(self, plane):
        if isinstance(self._observation, DerivedObservation):
            values = cc.find_keywords_in_headers(self._headers[1:], ['IMCMB'])
            repaired_values = _remove_processing_detritus(values, self._observation.observation_id)
            cc.update_plane_provenance_from_values(
                plane,
                self._repair_provenance_value,
                repaired_values,
                self._storage_name.collection,
                self._observation.observation_id,
            )

        if (
            ofr.is_processed(self._storage_name.file_name) or isinstance(self._observation, DerivedObservation)
        ) and 'jpg' not in self._storage_name.file_name:
            # not the preview artifact
            if plane.provenance is not None:
                plane.provenance.reference = (
                    f'http://archive.gemini.edu/searchform/filepre={self._storage_name.file_name}'
                )

    def _update_artifact(self, artifact):
        self._logger.debug(f'Begin _update_artifact for {artifact.uri}')
        for part in artifact.parts.values():
            for chunk in part.chunks:
                # JJK 20-01-25 - no cutouts
                chunk.naxis = None
                if chunk.time:
                    chunk.time_axis = None
                if chunk.position:
                    chunk.position_axis_1 = None
                    chunk.position_axis_2 = None
        self._logger.debug('End _update_artifact')


class MAROONXSpectralTemporal(MAROONXTemporal):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_energy_axis(1)
        # https://gemini.edu/instrumentation/maroon-x/capability
        # A dichroic splits the light into a
        # ‘blue’ (full range: 491-670nm, order 124-92;
        #       useful range: 499-663nm, order 122-93)
        # and
        # ‘red’ (649-920nm, order 94-67).
        for extension in range(0, len(self._headers)):
            bp.set('Chunk.energy.specsys', 'TOPOCENT', extension)
            camera_arm = self._headers[extension].get('ARM')
            if camera_arm and camera_arm in ['RED', 'BLUE']:
                bp.set('Chunk.energy.axis.axis.ctype', 'WAVE', extension)
                bp.set('Chunk.energy.axis.axis.cunit', 'nm', extension)
                bp.set('Chunk.energy.axis.range.start.pix', 0.5, extension)
                bp.set('Chunk.energy.axis.range.start.val', '_get_energy_chunk_range_start_val()', extension)
                bp.set('Chunk.energy.axis.range.end.pix', 1.5, extension)
                bp.set('Chunk.energy.axis.range.end.val', '_get_energy_chunk_range_end_val()', extension)
        # WF 07-01-25 resolving power ~85k
        bp.set('Chunk.energy.resolvingPower', 85000.0)
        self._logger.debug(f'End accumulate_blueprint')

    def _get_energy_chunk_range_end_val(self, ext):
        result = None
        if (camera_arm := self._headers[ext].get('ARM')) == 'RED':
            result = 920.0  # nm
        elif camera_arm == 'BLUE':
            result = 663.0  # nm
        else:
            raise mc.CadcException(f'Unexpected ARM value {camera_arm}.')
        return result

    def _get_energy_chunk_range_start_val(self, ext):
        result = None
        if (camera_arm := self._headers[ext].get('ARM')) == 'RED':
            result = 649.0  # nm
        elif camera_arm == 'BLUE':
            result = 499.0  # nm
        return result


class MAROONXSpatialSpectralTemporal(MAROONXSpectralTemporal):
    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((2, 3))
        for extension in range(0, len(self._headers)):
            self._accumulate_chunk_position_axes_blueprint(bp, extension)
        self._logger.debug(f'End accumulate_blueprint')


class Michelle(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')
        bp.configure_position_axes((1, 2))

    def _get_filter_name(self, ext):
        # DB - 08-04-19
        # Use header FILTERA and FILTERB values to determine the filter
        # bandpass
        # 08-04-19
        # Note that FILTERB sometimes has a value “Clear_B” and that should
        # be ignored.  Not sure if “Clear_A” is possible but maybe best to
        # allow for it.
        return self._search_through_keys(0, ['FILTER'])

    def _make_axes_consistent(self, chunk):
        super()._make_axes_consistent(chunk)
        if chunk.naxis == 4:
            if chunk.position is None:
                chunk.naxis = None
                if chunk.time_axis is not None:
                    chunk.time_axis = None
            else:
                chunk.naxis = 2
                chunk.time_axis = None

    def _reset_energy(self, observation_type, filter_name):
        # 'No Value' in filter_name test obs GN-2005A-C-14-45-002
        # DB 09-04-19 - “Blank-B” -> no energy.  It is a bias exposure
        # apparently and those shouldn't have energy WCS.
        # DB 23-04-19 - “No Value” -> no energy
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and ('Blank' in filter_name or 'No Value' in filter_name)
        )

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')

        # DB - 01-03-19
        #
        # no resolution for imaging mode
        #
        # Michelle spectroscopy:
        # json ‘disperser’, e.g. LowN, MedN1, MedN2, LowQ, Echelle
        # json ‘focal_plane_mask’, e.g. 2_pixels (need the ‘2’)
        #
        # Disp.    R    final R as function of slit width
        # --------------------------------------------------
        # LowN    200    R x 2/slit width
        # LowQ    110    R x 3/slit width
        # MedN1    1000    R x 2/slit width
        # MedN2    3000    R x 2/slit width
        # Echelle    30000    R x 2/slit width (very approximate)

        # DB - 04-04-19
        # The only solution for the michelle datasets (another ‘visitor’
        # instrument) would be to use filter info in the gsa..gsa_filters
        # table (as you do for PHOENIX).  No info to pass on to SVO folks to
        # add more filters.   Code would have to handle this case of two
        # filters with overlapping bandpasses.

        # DB - 15-04-19
        # Michelle:  From a hidden page of Michelle filters,
        # http://www.gemini.edu/sciops/instruments/michelle/imaging/filters,
        # hard-code QBlock and NBlock filters using values in that table
        # (ignoring the greater-than symbols)?  i.e. NBlock has central
        # bandpass of  (14+6.8)/2 microns and bandpass of 14-6.8 microns.
        # Ditto for QBlock:  (16.1+25)/2 and 25-16.1.  Grid_T should be
        # ignored for bandpass calculations but it would be good to keep it
        # in the filter name so in this case it would be F125B9 + Grid_T

        # 0 = R
        # 1 = ratio for slit width
        lookup = {
            'LowN': [200, 2.0],
            'LowQ': [110, 3.0],
            'MedN1': [1000, 2.0],
            'MedN2': [3000, 2.0],
            'Echelle': [30000, 2.0],
        }

        # 0 - central
        # 1 - lower
        # 2 - upper
        #
        # units are microns
        #
        # select filter_id,wavelength_central,
        # wavelength_lower,wavelength_upper from gsa..gsa_filters where
        # instrument="michelle";
        michelle = {
            'F103B10': [10.282990, 9.875740, 10.690240],
            'F105B53': [10.500000, 7.850000, 13.150000],
            'F112B21': [11.299610, 10.236910, 12.362310],
            'F116B9': [11.687970, 11.254020, 12.121920],
            'F125B9': [12.493580, 12.097480, 12.889680],
            'F128B2': [12.800000, 12.650000, 12.950000],
            'F14SA': [10.000000, 1.000000, 14.000000],
            'F161L': [20.500000, 16.100000, 25.000000],
            'F185B9B': [18.113860, 17.385410, 18.842310],
            'F198B27': [20.696080, 17.945430, 23.446730],
            'F209B42L': [20.900000, 16.500000, 25.300000],
            'F209B42S': [20.900000, 16.500000, 25.300000],
            'F66LA': [10.000000, 6.600000, 99.900000],
            'F66LB': [10.000000, 6.600000, 99.900000],
            'F79B10': [7.715890, 7.392540, 8.039240],
            'F88B10': [8.821500, 8.469400, 9.173600],
            'F97B10': [9.688450, 9.253750, 10.123150],
            'QBlock': [(16.1 + 25) / 2, 16.1, 25.0],
            'NBlock': [(14 + 6.8) / 2, 6.8, 14.0],
        }

        self.fm = self._filter_cache.get_filter_metadata(
            self._instrument,
            filter_name,
            self._lookup.telescope(self._storage_name.file_uri),
        )
        if self.fm is None:  # means filter_name not found
            # DB 09-04-19 - Use 100 microns for the initial max for michelle.
            w_max, w_min = self._multiple_filter_lookup(filter_name, michelle, wl_max=100)
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(w_max, w_min)
            self.fm.set_central_wl(w_max, w_min)
        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'Spectral WCS spectrum for {self._storage_name.obs_id}.')
            if data_product_type == DataProductType.SPECTRUM:
                disperser = self._lookup.disperser(self._storage_name.file_uri)
                focal_plane_mask = self._lookup.focal_plane_mask(self._storage_name.file_uri)
                slit_width = float(focal_plane_mask.split('_')[0])
                if disperser not in lookup:
                    self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                    raise mc.CadcException(
                        f'{self._instrument}: Mystery disperser {disperser} for {self._storage_name.obs_id}'
                    )
                self.fm.resolving_power = lookup[disperser][0] * lookup[disperser][1] / slit_width
        elif data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'Spectral WCS imaging mode for {self._storage_name.obs_id}.')
        else:
            raise mc.CadcException(
                f'{self._instrument}: no Spectral WCS support when '
                f'DataProductType {data_product_type} for '
                f'{self._storage_name.obs_id}'
            )

        # use the json value for bandpass_name value - it's representative of
        # multiple filters
        filter_name = self._lookup.filter_name(self._storage_name.file_uri)
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')


class Nici(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))

    def _get_filter_name(self, ext):
        # NICI - use filter names from headers, because there's a different
        # filter/header, and the JSON summary value obfuscates that
        return self._search_through_keys(ext, ['FILTER'])

    def _reset_energy(self, observation_type, filter_name):
        # DB 04-04-19
        # If one of the NICI filters is ‘Block’ then energy WCS should be
        # ignored for that extension.
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and 'Block' in filter_name
        )

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        # DB 22-08-19
        # Duplicate the line “‘FeII’: [1.644000, 1.631670, 1.656330]” in the
        # NICI filters replacing “FeII” with “[FeII]“?   Although it would be
        # better to have only [FeII] show up in the pick list and not both
        # Fe II and [Fe II]…

        # select filter_id, wavelength_central, wavelength_lower,
        #        wavelength_upper
        # from gsa..gsa_filters where instrument = 'NICI'
        # 0 - central
        # 1 - lower
        # 2 - upper
        #
        # dict with the barcodes stripped from the names as returned by query
        # DB - 22-02-19 - units are microns
        nici_lookup = {
            'Br-gamma': [2.168600, 2.153900, 2.183300],
            'CH4-H1%L': [1.628000, 1.619300, 1.636700],
            'CH4-H1%S': [1.587000, 1.579500, 1.594500],
            'CH4-H1%Sp': [1.603000, 1.594900, 1.611100],
            'CH4-H4%L': [1.652000, 1.619000, 1.685000],
            'CH4-H4%S': [1.578000, 1.547000, 1.609000],
            'CH4-H6.5%L': [1.701000, 1.652400, 1.749600],
            'CH4-H6.5%S': [1.596000, 1.537250, 1.654750],
            'CH4-K5%L': [2.241000, 2.187500, 2.294500],
            'CH4-K5%S': [2.080000, 2.027500, 2.132500],
            'FeII': [1.644000, 1.631670, 1.656330],
            'H2-1-0-S1': [2.123900, 2.110800, 2.137000],
            'H20-Ice-L': [3.090000, 3.020000, 3.150000],
            'H': [1.650000, 1.490000, 1.780000],
            'J': [1.250000, 1.150000, 1.330000],
            'K': [2.200000, 2.030000, 2.360000],
            'Kcont': [2.271800, 2.254194, 2.289406],
            'Kprime': [2.120000, 1.950000, 2.300000],
            'Ks': [2.150000, 1.990000, 2.300000],
            'Lprime': [3.780000, 3.430000, 4.130000],
            'Mprime': [4.680000, 4.550000, 4.790000],
        }
        filter_name = filter_name.split('_G')[0]
        filter_name = filter_name.replace('[FeII]', 'FeII')
        self.fm = self._filter_cache.get_filter_metadata(
            self._instrument,
            filter_name,
            self._lookup.telescope(self._storage_name.file_uri),
        )
        if data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging mode for {self._storage_name.obs_id}.')
            if self.fm is None:  # means filter_name not found
                w_max, w_min = self._multiple_filter_lookup(filter_name, nici_lookup)
                self.fm = svofps.FilterMetadata()
                self.fm.set_bandpass(w_max, w_min)
                self.fm.set_central_wl(w_max, w_min)
                self.fm.set_resolving_power(w_max, w_min)

            temp = self._lookup.filter_name(self._storage_name.file_uri)
            # NICI has two different bandpass names (most of the time) in
            # its two chunks.  Pat says in this case nothing will be put in
            # the bandpass name for the plane.  Add code to combine the two
            # chunk bandpass names to create a plane bandpass name only for
            # this instrument
            filter_name = temp.replace('[FeII]', 'FeII')
            self._build_chunk_energy(chunk, filter_name)
        else:
            raise mc.CadcException(
                f'{self._instrument}: Do not understand DataProductType '
                f'{data_product_type} from {self._storage_name.obs_id}'
            )
        self._logger.debug('End _update_energy')


class Nifs(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')

    def get_cd11(self, ext):
        # DB - 05-03-19 - NIFS needs a division by NAXIS1/2 for the
        # cdelta1/2 calculations.
        return RADIUS_LOOKUP[self._instrument] / self._headers[ext].get('NAXIS1')

    def get_cd22(self, ext):
        return self.get_cd11(ext)

    def _reset_energy(self, observation_type, filter_name):
        # filter_name == '', test obs is GN-CAL20050301-17-001
        # DB 18-04-19
        #
        # Gemini archive shows WaveBand=INVALID - ignore energy WCS
        # DB 23-04-19 - ‘no value’ -> no energy

        # DB 23-04-19 - Looks like NIRI observation GN-CAL20020623-1-011 is
        # skipping the ‘invalid’ values in FILTER1/2 headers.  So energy
        # should likely be skipped.  In this case the json value is
        # INVALID&INVALID.
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and ('Blocked' in filter_name or 'INVALID' in filter_name)
        )

    def _reset_position(self, observation_type, artifact_type):
        result = super()._reset_position(observation_type, artifact_type)

        # DB - 08-04-19 - json ra/dec values are null for
        # the file with things set to -9999.  Ignore
        # spatial WCS for these cases.

        # get the values from JSON directly, because the
        # function uses header values, which are set to
        # unlikely defaults
        ra = self._lookup.ra(self._storage_name.file_uri)
        dec = self._lookup.dec(self._storage_name.file_uri)
        if ra is None and dec is None:
            result = True
        elif ra is not None and math.isclose(ra, 0.0) and dec is not None and math.isclose(dec, 0.0):
            result = True
        return result

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')

        # Imaging: treated as for other instruments.
        #
        # Spectroscopy:
        #
        # Need grating name from json ‘disperser’ value
        # Need central wavelength from json ‘central_wavelength’ value
        # Need filter name from json ‘filter’ value
        # Need header NAXIS1 value in extension for number of pixels in
        # dispersion direction
        #
        # Then the top table on this page to create a lookup for
        # upper/lower wavelength values and spectral resolution for
        # the given grating/filter combination:
        # https://www.gemini.edu/sciops/instruments/nifs/ifu-spectroscopy/gratings.
        #
        # cdelta will be as for other instruments (upper - lower)/naxis

        # Use the filter names WITHOUT the ‘_298K’ suffix (since those are for
        # the cold filters at operating temperature).

        # key = grating name
        # key = associated filter
        # 0 = Central Wavelength (microns),
        # 1, 2 = Spectral Range,
        # 3 = Spectral Resolution,
        # 4 = Velocity Resolution (km/s)

        # From the second table:
        # Table 2: The NIFS gratings can be tuned to different central
        # wavelengths. The short and long limits of the possible tuned
        # central wavelengths and the associated filters required are:
        # Grating, Name, Short Wavelength Limit (microns), Long Wavelength Limit (microns)	Short Wavelength Filter	Long Wavelength Filter
        # Z	0.94	1.16	ZJ	ZJ
        # J	1.14	1.36	ZJ	JH
        # H	1.48	1.82	JH	HK
        # K	1.98	2.41	HK	HK

        # DB - 01-03-19
        # NIFS:  this apparently uses ‘K_Short’ and ‘K_Long’ gratings (in json
        # ‘disperser’) that I think are the same as the other K grating but tuned
        # (rotated) to different wavelengths than normal.  So these should use
        # the same resolution values as the K grating lookup.

        # DB - 04-04-19
        # NIFS gratings can be tuned (rotated) for slightly different
        # wavelengths.  But these lines should be added to the NIFS lookup:
        #
        # ‘J’: {‘JH’: [1.25, 1.15, 1.33, 6040.0, 49.6]}
        # ‘H’: {‘HK’: [1.65, 1.49, 1.80, 5290.0, 56.8]}
        #
        # i.e. the same values as the other J/H grating entries except for
        # the filter names.  What’s important are the 2nd and 3rd numbers
        # that are used to determine the bandpass since you use the json
        # ‘central_wavelength’ to establish that value.  That central
        # wavelength combined with the bandpass should set different
        # upper/lower wavelength limits for NIFS observations with the same
        # grating/filter combination.
        #
        # For the NIFS JH filter observations with disperser K those appear
        # to be observations that they use to measure noise or dark signal or
        # perhaps point to a very bright star for an acquisition observation.
        # There are other combinations I see, e.g. ZJ with K.  Those should
        # skip energy WCS calculation.  i.e. any filter/disperser combination
        # not in the lookup.   Also, I see a ‘Blocked’ filter for NIFS.
        # Those should also skip energy WCS.

        nifs_lookup = {
            'Z': {
                'ZJ': [1.05, 0.94, 1.15, 4990.0, 60.1],
            },
            'J': {
                'ZJ': [1.25, 1.15, 1.33, 6040.0, 49.6],
                'JH': [1.25, 1.15, 1.33, 6040.0, 49.6],
            },
            'H': {
                'JH': [1.65, 1.49, 1.80, 5290.0, 56.8],
                'HK': [1.65, 1.49, 1.80, 5290.0, 56.8],
            },
            'K': {
                'HK': [2.20, 1.99, 2.40, 5290.0, 56.7],
            },
            'K_Short': {
                'HK': [2.20, 1.98, 2.41, 5290.0, 56.7],
            },
            'K_Long': {
                'HK': [2.20, 1.98, 2.41, 5290.0, 56.7],
            },
        }

        self.fm = self._filter_cache.get_filter_metadata(
            self._instrument,
            filter_name,
            self._lookup.telescope(self._storage_name.file_uri),
        )

        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'spectroscopy for {self._storage_name.obs_id}.')
            grating = self._lookup.disperser(self._storage_name.file_uri)
            if grating in nifs_lookup:
                if filter_name in nifs_lookup[grating]:
                    self.fm = svofps.FilterMetadata('NIFS')
                    self.fm.set_bandpass(
                        nifs_lookup[grating][filter_name][2],
                        nifs_lookup[grating][filter_name][1],
                    )
                    self.fm.central_wl = self._lookup.central_wavelength(self._storage_name.file_uri)
                    self.fm.resolving_power = nifs_lookup[grating][filter_name][3]
                    if self.fm.central_wl is None:
                        # DB 14-04-21
                        #  if no central wavelength in json then use filter
                        #  bounds
                        self.fm.central_wl = nifs_lookup[grating][filter_name][0]
                        self._logger.warning(
                            f'JSON central_wavelength is None for '
                            f'{self._storage_name.obs_id}. Using '
                            f'{self.fm.central_wl} instead.'
                        )
                else:
                    self._logger.info(
                        f'No energy. filter_name {filter_name} with '
                        f'disperser {grating} for {self._storage_name.obs_id}'
                    )
                    fm = None
            else:
                self._logger.info(f'No energy. grating {grating} for {self._storage_name.obs_id}')
                fm = None
        elif data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'NIFS: imaging for {self._storage_name.obs_id}.')
            # DB - 01-03-19
            # NIFS images should just use the standard imaging procedure
            # for resolution (central_wavelength/bandpass).
        else:
            raise mc.CadcException(
                f'{self._instrument}: DataProductType {data_product_type} for {self._storage_name.obs_id}'
            )
        if self.fm is None or self.fm.central_wl is None:
            cc.reset_energy(chunk)
        else:
            self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End update_energy')

    def _update_position(self, part, chunk, ext):
        # DB - 01-18-19 - NIFS has no WCS info in
        # extension; use primary header
        #
        # DB - 04-03-19 - NIFS spatial WCS info in the
        # header has way too large a FOV so hard-code
        # this to the instrument's tiny 3" x 3" FOV.

        self._logger.debug('Begin update_position')
        n_axis1 = self._headers[-1]['NAXIS1']
        n_axis2 = self._headers[-1]['NAXIS2']

        # DB - 18-02-19 - for hard-coded field of views use:
        # CRVAL1  = RA value from json or header (degrees
        # CRVAL2  = Dec value from json or header (degrees)
        # CDELT1  = 5.0/3600.0 (Plate scale along axis1 in degrees/pixel
        #           for 5" size)
        # CDELT2  = 5.0/3600.0
        # CROTA1  = 0.0 / Rotation in degrees
        # NAXIS1 = 1
        # NAXIS2 = 1
        # CRPIX1 = 1.0
        # CRPIX2 = 1.0
        # CTYPE1 = RA---TAN
        # CTYPE2 = DEC--TAN

        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(0)
        header['CRVAL2'] = self.get_dec(0)
        header['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        header['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        header['CROTA1'] = 0.0
        # DB 05-03-19 - persist NAXIS values for NIFS
        header['NAXIS1'] = n_axis1
        header['NAXIS2'] = n_axis2
        # DB 18-05-21
        # EQUINOX and RADESYS values should be fine for all NIFS files since
        # that’s what is used for all other Gemini data as far as I’m aware.
        header['EQUINOX'] = 2000.0
        header['RADESYS'] = 'FK5'

        header['CRPIX1'] = self.get_crpix1(0)
        header['CRPIX2'] = self.get_crpix2(0)
        header['CD1_1'] = self.get_cd11(0)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(0)

        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, ext)
        if chunk is None:
            chunk = Chunk()
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        self._logger.debug('End update_position')


class Niri(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))

    def _get_filter_name(self, ext):
        # NIRI - prefer header keywords
        return self._search_through_keys(0, ['FILTER'])

    def _reset_energy(self, observation_type, filter_name):
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and (filter_name == '' or 'INVALID' in filter_name)
        )

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        reset_energy = False

        # values from
        # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/grisms
        # units are? page is in microns

        # DB 21-05-19
        # For lack of info on the slit for these odd observations use the
        # same resolution for the 4-pixel slit:  i.e. add “‘f6-cam’: 610.0”
        # for the J grism lookup (of which there may only be two cases)
        # and “‘f6-cam’: 780.0" for the K grism.  I don’t see any similar
        # types of observations for other grisms that are not ‘blank’
        # observations and so don’t have energy WCS.

        # DB 22-08-19
        # Guesstimate from changes in values for L filter/grism centered and
        # blue combinations.
        # M: 1100 for f6-2pixBl, 860 for f6-4pixBl and 490 for f6-6pixBl.
        # Skip energy for INVALID f-ratio.

        NIRI_RESOLVING_POWER = {
            'J': {
                'f6-2pix': 770.0,
                'f6-4pix': 610.0,
                'f6-cam': 610.0,
                'f6-6pix': 460.0,
                'f6-2pixBl': 770.0,
                'f6-4pixBl': 650.0,
                'f6-6pixBl': 480.0,
                'f32-4pix': 1000.0,
                'f32-6pix': 620.0,  # f32-7pix
                'f32-9pix': 450.0,  # f32-10pix
            },
            'H': {
                'f6-2pix': 1650.0,
                'f6-4pix': 825.0,
                'f6-6pix': 520.0,
                'f6-2pixBl': 1650.0,
                'f6-4pixBl': 940.0,
                'f6-6pixBl': 550.0,
                'f32-4pix': 880.0,
                'f32-6pix': 630.0,  # f32-7pix
                'f32-9pix': 500.0,  # f32-10pix
            },
            'L': {
                'f6-2pix': 1100.0,
                'f6-4pix': 690.0,
                'f6-6pix': 460.0,
                'f6-2pixBl': 1100.0,
                'f6-4pixBl': 770.0,
                'f6-6pixBl': 490.0,
            },
            'M': {
                'f6-2pix': 1100.0,
                'f6-4pix': 770.0,
                'f6-6pix': 460.0,
                'f6-2pixBl': 1100.0,
                'f6-4pixBl': 860.0,
                'f6-6pixBl': 490.0,
            },
            'K': {
                'f6-2pix': 1300.0,
                'f6-4pix': 780.0,
                'f6-cam': 780.0,
                'f6-6pix': 520.0,
                'f6-2pixBl': 1300.0,
                'f6-4pixBl': 780.0,
                'f6-6pixBl': 520.0,
                'f32-4pix': 1280.0,
                'f32-6pix': 775.0,  # f32-7pix
                'f32-9pix': 570.0,  # f32-10pix
            },
        }

        # NIRI.Hgrismf32-G5228

        # No energy information is determined for darks.  The
        # latter are sometimes only identified by a 'blank' filter.  e.g.
        # NIRI 'flats' are sometimes obtained with the filter wheel blocked
        # off.

        # The focal_plane_mask has values like f6-2pixBl (2-pixel wide blue
        # slit used with f/6 camera) or f32-6pix (6-pixel wide slit [centered]
        # with f/32 camera).   So f#-#pix with or without ‘Bl’.   But there
        # seems to be some inconsistency in the above web page and actual
        # slits.  I don’t think there are 7- or 10-pixel slits used with the
        # f32 camera.  The Gemini archive pull-down menu only gives f32-6pix
        # and f32-9pix choices.  Assume these refer to the 7 and 10 pixel
        # slits in the web page. Use the ‘disperser’ value to lookup the
        # ‘Estimated Resolving Power’ for the slit in the beam given by the
        # focal_plane_mask value.

        # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/blocking-filters

        # DB - 08-09-19 - NIRI data with ‘pinhole’:  ignore energy WCS for
        # these for any spectroscopic observations.  They are used for image
        # quality tests and instrument alignment.  For imaging observations
        # these shouldn’t impact anything (I hope).  e.g. does this one work?
        # GN-2005A-DD-12-12-017

        # DB - 15-04-19
        # For r_ratio:
        # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/grisms,
        # parenthetically states “w. either slits” for the f/6 K section.  So
        # duplicate lines 1268-1270 of main_app.py with “Bl” appended after
        # ‘pix’.

        if 'Jcon(112)_G0235' in filter_name:
            # DB - 01-04-19 The G0235 filter is listed as ‘damaged’ on the
            # Gemini NIRI filters web site:
            # https://www.gemini.edu/sciops/instruments/niri/imaging/filters.
            # Not enough info is given there for SVO to add this filter to
            # their system.  Hardcode a central wavelength of 1.1232 microns
            # and a FWHM of 0.0092 microns
            filter_md = svofps.FilterMetadata('NIRI')
            filter_md.central_wl = 1.1232
            filter_md.bandpass = 0.0092
        elif 'Msort' in filter_name or 'Mgrism' in filter_name:
            # DB - 15-04-19
            # The NIRI blocking filter page,
            # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/blocking-
            # filters, gives no transmission data for this filter for the SVO
            # to use.  Will need to hard-code it to central wavelength of
            # (4.4+6)/2 and width of 6-4.4 microns.

            filter_md = svofps.FilterMetadata('NIRI')
            filter_md.set_central_wl(6.0, 4.4)
            filter_md.set_bandpass(6.0, 4.4)
        else:
            filter_md = self._filter_cache.get_filter_metadata(
                self._instrument,
                filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )
            if filter_md is None:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(
                    f'{self._instrument}: mystery filter {filter_name} for {self._storage_name.obs_id}'
                )

        filter_name = self._lookup.filter_name(self._storage_name.file_uri)
        if data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging for {self._storage_name.obs_id}.')
            self.fm = filter_md
            self.fm.adjust_resolving_power()
        elif data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS spectroscopy for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata('NIRI')
            self.fm.central_wl = filter_md.central_wl
            self.fm.bandpass = filter_md.bandpass
            # add the 'split' call because NIRI: Mystery disperser value
            # Mgrism_G5206 for GN-2007B-Q-75-61-003
            disperser = self._lookup.disperser(self._storage_name.file_uri).split('_')[0]
            if disperser in [
                'Jgrism',
                'Jgrismf32',
                'Hgrism',
                'Hgrismf32',
                'Kgrism',
                'Kgrismf32',
                'Lgrism',
                'Mgrism',
            ]:
                bandpass_name = disperser[0]
                f_ratio = self._lookup.focal_plane_mask(self._storage_name.file_uri)
                self._logger.debug(
                    f'Bandpass name is {bandpass_name} f_ratio is {f_ratio} for {self._storage_name.obs_id}'
                )
                if bandpass_name in NIRI_RESOLVING_POWER and f_ratio in NIRI_RESOLVING_POWER[bandpass_name]:
                    self.fm.resolving_power = NIRI_RESOLVING_POWER[bandpass_name][f_ratio]
                elif 'pinhole' in f_ratio:
                    self._logger.info(f'Pinhole. Setting energy to None for {self._storage_name.obs_id}')
                    reset_energy = True
                elif f_ratio == 'INVALID':
                    self._logger.info(f'INVALID f_ratio. Setting energy to None for {self._storage_name.obs_id}')
                    reset_energy = True
                else:
                    self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                    raise mc.CadcException(
                        f'{self._instrument} Mystery bandpass name '
                        f'{bandpass_name} or f_ratio {f_ratio} for '
                        f'{self._storage_name.obs_id}.'
                    )
            else:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(
                    f'{self._instrument} Mystery disperser value {disperser} for {self._storage_name.obs_id}'
                )
        else:
            raise mc.CadcException(
                f'{self._instrument}: Do not understand mode {data_product_type} for {self._storage_name.obs_id}'
            )

        if reset_energy:
            cc.reset_energy(chunk)
        else:
            self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_position(self, part, chunk, extension):
        self._logger.info(f'Begin _update_position')
        # DB 07-06-21
        # The extension CD values that are very, very close to 0 cause the
        # problems with Spatial WCS:
        # ERROR: spherepoly_from_array: a line segment overlaps or polygon too
        #        large
        # Try to use the primary values if this error occurs - there's an extra
        # '5' in the exponent
        if len(self._headers) > 1:
            pdu = self._headers[0]
            hdu0 = self._headers[1]
            pdu_cd1_1 = pdu.get('CD1_1')
            hdu0_cd1_1 = hdu0.get('CD1_1')
            if (not math.isclose(pdu_cd1_1, hdu0_cd1_1) and math.isclose(pdu_cd1_1 * 1e-50, hdu0_cd1_1)) or (
                hdu0.get('CRPIX1') != 0.0
                and hdu0.get('CRVAL1') != 0.0
                and hdu0.get('CRPIX2') != 0.0
                and hdu0.get('CRVAL2') != 0.0
            ):
                pdu['NAXIS1'] = hdu0.get('NAXIS1')
                pdu['NAXIS2'] = hdu0.get('NAXIS2')
                wcs_parser = FitsWcsParser(pdu, self._storage_name.obs_id, extension)
                if chunk is None:
                    chunk = Chunk()
                    part.chunks.append(chunk)
                wcs_parser.augment_position(chunk)
                if chunk.position is not None:
                    chunk.position_axis_1 = 1
                    chunk.position_axis_2 = 2
                    chunk.position.coordsys = pdu.get('FRAME')
                    chunk.position.equinox = mc.to_float(pdu.get('EQUINOX'))
            else:
                if chunk is not None and chunk.position is not None:
                    cc.reset_position(chunk)
        self._logger.info('End _update_position')


class Oscir(GeminiMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')

    def get_art_product_type(self, ext):
        # TODO - don't know what else to do right now
        return ProductType.SCIENCE

    def get_calibration_level(self, ext):
        result = CalibrationLevel.RAW_STANDARD
        if self._storage_name.file_id.startswith('r'):
            # DB 23-02-21
            # The best thing to do with OSCIR 'r' files is to add them as a
            # second cal level 2 plane and use the same metadata as the
            # unprocessed plane.
            result = CalibrationLevel.CALIBRATED
        return result

    def get_cd11(self, ext):
        return RADIUS_LOOKUP[self._instrument]

    def get_cd22(self, ext):
        return self.get_cd11(ext)

    def get_crpix1(self, ext, keyword='NAXIS1'):
        naxis1 = self._headers[ext].get(keyword)
        if naxis1 is None:
            result = None
        else:
            result = naxis1 / 2.0
        return result

    def get_crpix2(self, ext):
        return self.get_crpix1(ext, 'NAXIS2')

    def get_dec(self, ext):
        ra_ignore, dec = self.get_sky_coord(ext, 'RA_TEL', 'DEC_TEL')
        return dec

    def get_exposure(self, ext):
        # DB - 20-02-19 - json ‘exposure_time’ is in minutes, so multiply
        # by 60.0.
        return self._lookup.exposure_time(self._storage_name.file_uri) * 60.0

    def get_ra(self, ext):
        ra, dec_ignore = self.get_sky_coord(ext, 'RA_TEL', 'DEC_TEL')
        return ra

    def get_time_function_val(self, ext):
        # DB - 20-02-19 - OSCIR json ‘ut_datetime’ is not correct.  Must
        # concatenate DATE-OBS and UTC1 values and convert to MJD as usual
        # or use MJD directly (seems to be correct starting MJD)
        #
        # DB - 06-03-19 json ut_date_time Gemini is based on DATE keyword.
        # Better to use MJD header keyword value directly.
        return self._headers[ext].get('MJD')

    def _make_axes_consistent(self, chunk):
        # DB - 01-02-21
        # OSCIR data should all have NAXIS=6 in the header. Axis 1/2 are
        # position axes.  Axis 3 is the number of chop positions (1 or usually
        # = 2), axis 5 is the number of nod positions (1 or usually = 2), axis
        # 4 gives the number of ‘savesets’ per nod position and axis 6 gives
        # the number of ‘nod sets’.
        #
        # For that particular example you gave me I think the NAXIS6 = 1 is
        # incorrect.  Think it should equal the value of NODSETS.  And TOTFRMS
        # = 360 should equal naxis1 x naxis2 x naxis3 x naxis4.
        #
        # Spatial cutouts would use axes 1 and 2, although the precise
        # positions of each chop/nod position are not captured by the CAOM2
        # data.
        super()._make_axes_consistent(chunk)
        if chunk.naxis is not None and chunk.naxis == 6:
            chunk.naxis = 2
            chunk.time_axis = None

    def _update_energy(self, chunk, data_product_type, filter_name):
        """oscir-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin _update_energy')

        # Filter info here:
        # http://www.gemini.edu/sciops/instruments/oscir/oscirFilterList.html
        # No filter provided in json; use FILTER keyword.
        # e.g. ‘S_12.5 (-11775)’ = ‘12.5’ in table.
        # It looks like only the 'r' files have filter ids.

        # DB 23-04-19
        # Split on whitespace and name the filters S_8.8 and IHW_(17-19)
        # respectively.  And change the keys of the lookup table to these
        # values as well.Randomly looking at observations these are other
        # filter values:  S_12.5 (-11774), S_11.7 (-22125), S_7.9 (-63440),
        # and N_wide (-1154).   The latter instead of the lookup “N”.  So I’m
        # guessing the first six key values in oscir_lookup should be
        # S_7.9, S_8.8, S_9.8, S_11.7 and S_12.5.   Then N_wide, IHW_(17-19).
        #
        # 0 - central wavelenth
        # 1 - bandpass
        # units are microns
        oscir_lookup = {
            'S_7.9': [7.91, 0.755],
            'S_8.8': [8.81, 0.871],
            'S_9.8': [9.80, 0.952],
            'S_10.3': [10.27, 1.103],
            'S_11.7': [11.70, 1.110],
            'S_12.5': [12.49, 1.156],
            'N_wide': [10.75, 5.230],
            'IHW_(17-19)': [18.17, 1.651],
            'Q3': [20.8, 1.650],
        }

        temp = filter_name
        filter_name = temp.split()[0]
        if filter_name not in oscir_lookup:
            self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
            raise mc.CadcException(
                f'{self._instrument}: Mystery FILTER keyword {filter_name} for {self._storage_name.obs_id}'
            )
        # DB 23-02-21
        # OSCIR files with the 'r' prefix will cause issues because they are
        # classified as spectra so remove the check
        #   'if data_product_type == DataProductType.IMAGE:'
        # and assume they are all images.  The filter bandpass is the best we
        # can do for the energy WCS in any case.
        self._logger.debug(f'SpectralWCS imaging mode for {self._storage_name.obs_id}.')
        self.fm = svofps.FilterMetadata()
        self.fm.central_wl = oscir_lookup[filter_name][0]
        self.fm.bandpass = oscir_lookup[filter_name][1]
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug(f'End _update_energy')

    def _update_position(self, part, chunk, ext):
        header = self._headers[ext]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(ext)
        header['CRVAL2'] = self.get_dec(ext)
        header['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        header['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        header['CROTA1'] = 0.0
        header['CRPIX1'] = self.get_crpix1(ext)
        header['CRPIX2'] = self.get_crpix2(ext)
        header['CD1_1'] = self.get_cd11(ext)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(ext)

        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, ext)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2

        chunk.position.coordsys = header.get('FRAMEPA')
        chunk.position.equinox = mc.to_float(header.get('EQUINOX'))
        self._logger.debug('End _update_position')


class Phoenix(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)

    def get_art_product_type(self, ext):
        result = ProductType.SCIENCE
        object_value = self._headers[ext].get('OBJECT').lower()
        if (
            'flat ' in object_value
            or 'dark ' in object_value
            or 'arc' in object_value
            or 'comp' in object_value
            or 'lamp' in object_value
            or 'comparison' in object_value
            or 'slit' in object_value
        ):
            result = ProductType.CALIBRATION
        return result

    def get_calibration_level(self, ext):
        result = CalibrationLevel.RAW_STANDARD
        if self._storage_name.file_id.lower().startswith('p'):
            result = CalibrationLevel.CALIBRATED
        return result

    def get_obs_type(self, ext):
        # DB - 18-02-19 - make PHOENIX searches more useful by making
        # estimated guesses at observation type
        #
        # Looking at a few random nights it looks like a reasonable way to
        # determine if an observation is a FLAT is to look for ‘flat’ (any
        # case combination) in the json ‘object’ value or ‘gcal’ or ‘GCAL’.
        # But if ‘gcal’ AND ‘dark’ are in the ‘object’ string it’s a DARK.
        # (see below)

        # Ditto for an ARC if json ‘object’ contains the string ‘arc’

        # It’s a DARK obs type if the value of FITS header VIEW_POS contains
        # the string ‘dark’ OR if ‘dark’ is in json ‘object’ string.
        # (There appear to be cases where darks are taken with the
        # VIEW_POS = open.)  When it’s a dark do NOT generate WCS for energy
        # since the filter is often ‘open’ and energy info isn’t important for
        # DARK exposures.

        result = 'OBJECT'
        object_value = self._lookup.object(self._storage_name.file_uri).lower()
        view_pos = self._headers[ext].get('VIEW_POS')
        if 'flat' in object_value:
            result = 'FLAT'
        elif 'dark' in object_value or (view_pos is not None and 'dark' in view_pos):
            result = 'DARK'
        elif 'gcal' in object_value:
            result = 'FLAT'
        elif 'arc' in object_value:
            result = 'ARC'
        elif 'slit' in object_value:
            # DB 22-02-21
            # These are images that show the slit location so I think it’s best
            # to add a new OBSTYPE of SLIT (sort of like MASK for GMOS-N/S).
            result = 'SLIT'
        return result

    def _get_filter_name(self, ext):
        filter_name = super()._get_filter_name(ext)
        if filter_name is None or len(filter_name.strip()) == 0:
            filter_name = self._search_through_keys(ext, ['CVF_POS', 'FILT_POS'])
        return filter_name

    def _make_axes_consistent(self, chunk):
        super()._make_axes_consistent(chunk)
        if chunk.naxis == 4:
            if chunk.position is None:
                chunk.naxis = None
                if chunk.time_axis is not None:
                    chunk.time_axis = None
            else:
                chunk.naxis = 2
                chunk.time_axis = None

        # DB - 05-06-20
        # That’s a composite observation (but with no way of determining the 4
        # members from the header) and it is an extracted spectrum and
        # (despite some header info suggesting otherwise) has no wavelength
        # scale. LINEAR must be a reference to the fact that the spacing of
        # each pixel is constant.  Since the scale is simply in pixels….
        # CTYPE1 will refer to the energy axis.   Would likely make more sense
        # for a value of PIXEL.
        ctype = self._headers[0].get('CTYPE1')
        if ctype is not None and ctype in [
            'LINEAR',
            'PIXEL',
        ]:
            chunk.naxis = None
            chunk.position_axis_1 = None
            chunk.position_axis_2 = None
            chunk.time_axis = None
            chunk.energy_axis = None

    def _reset_energy(self, observation_type, filter_name):
        # DB 07-06-21
        # set energy to None
        return super()._reset_energy(observation_type, filter_name) or (
            filter_name is not None and 'invalid' in filter_name
        )

    def _reset_position(self, observation_type, artifact_type):
        result = super()._reset_position(observation_type, artifact_type)

        # DB 30-04-19
        # Looks like many relatively recent PHOENIX files have no RA/Dec
        # values in the header and so will have no spatial WCS.
        # Base this decision on json null values.  But looking at all
        # of the PHOENIX data from 2016 until 3 December 2017 it looks
        # like json ra/dec values are either null or 0.0 for all.
        # In both cases spatial WCS should be ignored.  (It will be
        # very difficult for users to find anything of interest in
        # these datasets other than searching by free-form target
        # names…)  PHOENIX returned as a visitor instrument in May 2016
        # after about 5 years away.

        ra = self._lookup.ra(self._storage_name.file_uri)
        dec = self._lookup.dec(self._storage_name.file_uri)
        if ra is None and dec is None:
            result = True
        elif ra is not None and math.isclose(ra, 0.0) and dec is not None and math.isclose(dec, 0.0):
            result = True
        return result

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')
        # DB - 12-02-19 - Phoenix should be the same as TReCS but uses
        # NAXIS2 for the length of the dispersion axis.

        # select filter_id, wavelength_central, wavelength_lower,
        #     wavelength_upper
        # from gsa..gsa_filters where instrument = 'PHOENIX'
        # 0 - central
        # 1 - lower
        # 2 - upper
        #
        # units are microns
        #
        # dict with the filter wheels stripped from the names as returned by
        # query
        #
        # DB 30-04-19
        # A new filter not previously available.  Add info from here,
        # https://www.noao.edu/kpno/phoenix/filters.html

        PHOENIX = {
            '2030': [4.929000, 4.808000, 5.050000],
            '2150': [4.658500, 4.566000, 4.751000],
            '2462': [4.078500, 4.008000, 4.149000],
            '2734': [3.670500, 3.610000, 3.731000],
            '2870': [3.490500, 3.436000, 3.545000],
            '3010': [3.334500, 3.279000, 3.390000],
            '3100': [3.240000, 3.180000, 3.300000],
            '3290': [3.032500, 2.980000, 3.085000],
            '4220': [2.370000, 2.348000, 2.392000],
            '4308': [2.322500, 2.296000, 2.349000],
            '4396': [2.272500, 2.249000, 2.296000],
            '4484': [2.230000, 2.210000, 2.250000],
            '4578': [2.185000, 2.160000, 2.210000],
            '4667': [2.143000, 2.120000, 2.166000],
            '4748': [2.104000, 2.082000, 2.126000],
            '6073': [1.647000, 1.632000, 1.662000],
            '6420': [1.557500, 1.547000, 1.568000],
            '7799': [1.280500, 1.271000, 1.290000],
            '8265': [1.204500, 1.196000, 1.213000],
            '9232': [1.083000, 1.077000, 1.089000],
            'L2870': [3.490500, 3.436000, 3.545000],
            '9440': [1.058500, 1.053000, 1.064000],
            'open': [3.0, 1.0, 5.0],
        }

        # DB - 12-02-19 - Note that the parenthetical numbers
        # after the Phoenix filter names (in the header) indicates the
        # filter wheel slot the filter is in and may occasionally change
        # so should be disregarded.
        if len(filter_name) > 0:
            filter_name = filter_name.split()[0]
            # found some files with '_' in the name
            if '_' in filter_name:
                filter_name = filter_name.split('_')[0]

        self._logger.debug(f'filter_name is {filter_name} for {self._storage_name.obs_id}')
        self.fm = svofps.FilterMetadata('Phoenix')
        if data_product_type in [
            DataProductType.SPECTRUM,
            DataProductType.IMAGE,
        ]:
            self._logger.debug(f'DataProductType {data_product_type} for {self._storage_name.obs_id}.')
            if filter_name in PHOENIX:
                self.fm.set_bandpass(PHOENIX[filter_name][2], PHOENIX[filter_name][1])
                self.fm.central_wl = PHOENIX[filter_name][0]
            elif len(filter_name) == 0:
                # DB 11-02-21
                # With open filter in Phoenix the band pass coverage should
                # be from 1 to 5 microns, so central wavelength of 3 microns
                # and bandpass of 4 microns. Lines 2777 to 2778 should be
                # changed as well as adding an ‘open_(1)’ filter
                self.fm.set_bandpass(5.0, 1.0)
                self.fm.set_central_wl(5.0, 1.0)
            else:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(
                    f'{self._instrument} mystery filter name {filter_name} for {self._storage_name.obs_id}'
                )
        else:
            raise mc.CadcException(
                f'{self._instrument} Unsupported DataProductType '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End update_energy')

    def _update_position(self, part, chunk, ext):
        self._logger.debug('Begin _update_position')
        header = self._headers[ext]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(ext)
        header['CRVAL2'] = self.get_dec(ext)
        header['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        header['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = self.get_crpix1(ext)
        header['CRPIX2'] = self.get_crpix2(ext)
        header['CD1_1'] = self.get_cd11(ext)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(ext)
        temp = header.get('EQUINOX')
        if temp is None or math.isclose(temp, 0.0):
            header['EQUINOX'] = header.get('EPOCH')

        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, ext)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        self._logger.debug('End _update_position')


class Texes(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.clear('Plane.provenance.reference')
        bp.configure_position_axes((1, 2))

    def get_calibration_level(self, ext):
        result = CalibrationLevel.RAW_STANDARD
        if '_red' in self._storage_name.file_uri.lower() or '_sum' in self._storage_name.file_uri.lower():
            result = CalibrationLevel.CALIBRATED
        return result

    def get_cd11(self, ext, keyword='CDELT1'):
        return RADIUS_LOOKUP[self._instrument]

    def get_cd22(self, ext):
        return self.get_cd11(ext)

    def get_data_product_type(self, ext):
        return DataProductType.SPECTRUM

    def get_dec(self, ext):
        # bHROS, TEXES: ra/dec not in json
        return self._headers[ext].get('DEC')

    def get_exposure(self, ext):
        # DB - 07-03-19 -  Use header ‘OBSTIME’ value as exposure time in
        # seconds.
        return mc.to_float(self._headers[ext].get('OBSTIME'))

    def _get_filter_name(self, ext):
        return None

    def get_ra(self, ext):
        # bHROS, TEXES: ra/dec not in json
        return self._headers[ext].get('RA')

    def get_target_type(self, ext):
        return TargetType.OBJECT

    def _update_energy(self, chunk, data_product_type, filter_name):
        # DB - 07-03-19
        # TEXES Spectroscopy
        #
        # Some special code will be needed for datalabels/planes.  There are
        # no datalabels in the FITS header.  json metadata (limited) must be
        # obtained with URL like
        # https://archive.gemini.edu/jsonsummary/canonical/
        #     filepre=TX20170321_flt.2507.fits.
        # Use TX20170321_flt.2507 as datalabel.  But NOTE:  *raw.2507.fits and
        # *red.2507.fits are two planes of the same observation. I’d suggest
        # we use ‘*raw*’ as the datalabel and ‘*red*’ or ‘*raw*’ as the
        # appropriate product ID’s for the science observations.  The ‘flt’
        # observations do not have a ‘red’ plane.  The json document contains
        # ‘filename’ if that’s helpful at all.  The ‘red’ files do not exist
        # for all ‘raw’ files.
        #
        #
        # Header OBSTYPE appears to be correct; not json obs_type.
        #
        # No previews are generated by Gemini
        #
        # Energy WCS:
        #
        # Central wavelength is given by 10,000/header(WAVENO0).  I have to do
        # some more investigation to see if we can determine wavelength
        # coverage (i.e. see if I can identify the echelle/echelon info from
        # the header - I don’t think so).  For now use 0.25 microns as the
        # fixed FWHM bandpass.
        self._logger.debug('Begin _update_energy')
        if data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS spectral mode for {self._storage_name.obs_id}.')
            self.fm = svofps.FilterMetadata('TEXES')
            self.fm.central_wl = 10000 / self._headers[0].get('WAVENO0')
            self.fm.bandpass = 0.25
        else:
            # data_type/observing mode is always spectroscopy
            raise mc.CadcException(
                f'{self._instrument}: mystery data product type '
                f'{data_product_type} for {self._storage_name.obs_id}'
            )
        self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_position(self, part, chunk, ext):
        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(0)
        header['CRVAL2'] = self.get_dec(0)
        header['CDELT1'] = RADIUS_LOOKUP[self._instrument]
        header['CDELT2'] = RADIUS_LOOKUP[self._instrument]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = self.get_crpix1(0)
        header['CRPIX2'] = self.get_crpix2(0)
        header['CD1_1'] = self.get_cd11(0)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(0)
        wcs_parser = FitsWcsParser(header, self._storage_name.obs_id, 0)
        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        self._logger.debug('End _update_chunk_position')


class Trecs(GeminiMapping):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)

    def get_obs_intent(self, ext):
        obs_class = self._lookup.observation_class(self._storage_name.file_uri)
        result = ObservationIntentType.CALIBRATION
        if obs_class is None:
            obs_type = self.get_obs_type(ext)
            if obs_type in CAL_VALUES:
                result = ObservationIntentType.CALIBRATION
            else:
                result = ObservationIntentType.SCIENCE
                data_label = self._headers[ext].get('DATALAB')
                if data_label is not None and '-CAL' in data_label:
                    result = ObservationIntentType.CALIBRATION
        elif 'science' in obs_class:
            result = ObservationIntentType.SCIENCE
        return result

    def _make_axes_consistent(self, chunk):
        super()._make_axes_consistent(chunk)
        if chunk.naxis == 4:
            if chunk.position is None:
                chunk.naxis = None
                if chunk.time_axis is not None:
                    chunk.time_axis = None
            else:
                chunk.naxis = 2
                chunk.time_axis = None

    def _reset_energy(self, observation_type, filter_name):
        # DB 23-04-19
        # GMOS GS-2005A-Q-26-12-001.  Lots of missing metadata, including
        # release date so no energy.  Ditto for TReCS GS-CAL20041206-6-007.
        return super()._reset_energy(observation_type, filter_name) or (filter_name is not None and filter_name == '')

    def _update_energy(self, chunk, data_product_type, filter_name):
        self._logger.debug('Begin _update_energy')

        # Look at the json ‘disperser’ value.  If = LowRes-20" then
        # resolving power = 80.  If LowRes-10" then resolving power = 100.
        # There was a high-res mode but perhaps never used.  Again,
        # naxis = NAXIS1 header value and other wcs info as above for NIRI.
        # But might take some string manipulation to match filter names with
        # SVO filters.

        if filter_name == 'Datum+Datum':
            # DB 22-08-19
            # I think it might be better to not provide wavelength info rather
            # than a very wide range.  That’s how some other TReCS data with
            # null filter names are treated.
            self._logger.info(f'Filter is {filter_name}, no Spectral WCS for {self._storage_name.obs_id}')
        else:
            if filter_name == 'Qone-17.8um':
                # DB 22-08-19
                # Qone-17.8um filter.  I can find no info about the bandpass
                # for that filter on the web but this info was in the old
                # gsa_filters table. Hardcode lower/upper bandpasses of 17.3
                # and 18.17 microns.
                w_min = 17.3
                w_max = 18.17
                self.fm = svofps.FilterMetadata()
                self.fm.set_bandpass(w_max, w_min)
                self.fm.set_central_wl(w_max, w_min)
                self.fm.set_resolving_power(w_max, w_min)
            else:
                self.fm = self._filter_cache.get_filter_metadata(
                    self._instrument,
                    filter_name,
                    self._lookup.telescope(self._storage_name.file_uri),
                )
            if self.fm is None:
                self._observable.rejected.record(mc.Rejected.MYSTERY_VALUE, self._storage_name.file_name)
                raise mc.CadcException(f'{self._instrument}: Mystery filter {filter_name}')
            if data_product_type == DataProductType.IMAGE:
                self._logger.debug(f'imaging mode for {self._storage_name.obs_id}.')
            elif data_product_type == DataProductType.SPECTRUM:
                self._logger.debug(f'LS|Spectroscopy mode for {self._storage_name.obs_id}.')
                disperser = self._lookup.disperser(self._storage_name.file_uri)
                if disperser is not None:
                    if disperser == 'LowRes-20':
                        self.fm.resolving_power = 80.0
                    elif disperser == 'LowRes-10':
                        self.fm.resolving_power = 100.0
            else:
                raise mc.CadcException(
                    f'{self._instrument}: Do not understand mode '
                    f'{data_product_type} for {self._storage_name.obs_id}'
                )
            self._build_chunk_energy(chunk, filter_name)
        self._logger.debug('End _update_energy')

    def _update_position(self, part, chunk, ext):
        # DB 22-08-19
        # For the file with CTYPE1 = 0 in HDU 1, all
        # of the other WCS info appears to be fine
        # (i.e. identical with the primary header).
        # Use the HDU 0 values.
        if len(self._headers) > 1:
            ctype1 = self._headers[ext].get('CTYPE1')
            if isinstance(ctype1, str):
                # value repair for a small subset of TReCS files  :(
                # test is rS20060306S0090, GS-2005B-Q-10-63-003
                if ctype1 == '0':
                    self._headers[ext]['CTYPE1'] = 'RA---TAN'
                wcs_parser = FitsWcsParser(self._headers[ext], self._storage_name.obs_id, ext)
            else:
                wcs_parser = FitsWcsParser(self._headers[0], self._storage_name.obs_id, ext)
        else:
            wcs_parser = FitsWcsParser(self._headers[0], self._storage_name.obs_id, ext)

        if chunk is None:
            chunk = Chunk()
            part.chunks.append(chunk)

        wcs_parser.augment_position(chunk)
        if chunk.position is not None:
            chunk.position_axis_1 = 1
            chunk.position_axis_2 = 2
            chunk.position.coordsys = self._headers[0].get('FRAME')
            chunk.position.equinox = mc.to_float(self._headers[0].get('EQUINOX'))


def get_obs_type(lookup, storage_name):
    result = lookup.observation_type(storage_name.file_uri)
    obs_class = lookup.observation_class(storage_name.file_uri)
    if obs_class is not None and (obs_class == 'acq' or obs_class == 'acqCal'):
        result = 'ACQUISITION'
    return result


def mapping_factory(storage_name, clients, reporter, observation, config):
    metadata_lookup = GeminiMetadataLookupStorageName(storage_name)
    provenance_finder = ProvenanceFinder(clients, config)
    inst = metadata_lookup.instrument(storage_name.file_uri)
    lookup = {
        Inst.BHROS: Bhros,
        Inst.CIRPASS: Cirpass,
        Inst.F2: F2,
        Inst.FLAMINGOS: Flamingos,
        Inst.ALOPEKE: Fox,
        Inst.ZORRO: Fox,
        Inst.GMOS: Gmos,
        Inst.GMOSS: Gmos,
        Inst.GMOSN: Gmos,
        Inst.GNIRS: Gnirs,
        Inst.GPI: Gpi,
        Inst.GSAOI: Gsaoi,
        Inst.HOKUPAA: Hokupaa,
        Inst.HRWFS: Hrwfs,
        Inst.IGRINS: Igrins,
        Inst.IGRINS2: Igrins2,
        Inst.MICHELLE: Michelle,
        Inst.NICI: Nici,
        Inst.NIFS: Nifs,
        Inst.NIRI: Niri,
        Inst.OSCIR: Oscir,
        Inst.PHOENIX: Phoenix,
        Inst.TEXES: Texes,
        Inst.TRECS: Trecs,
    }
    result = None
    if inst in lookup:
        result = lookup.get(inst)(
            storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
        )
    elif inst is Inst.GHOST:
        obs_type = get_obs_type(metadata_lookup, storage_name)
        if obs_type in ['ARC', 'BIAS', 'FLAT']:
            result = GHOSTSpectralTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
        else:
            result = GHOSTSpatialSpectralTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
    elif inst is Inst.GRACES:
        obs_type = get_obs_type(metadata_lookup, storage_name)
        mode = metadata_lookup.mode(storage_name.file_uri)
        # DB 30-04-19
        # Ignore spatial WCS for any GRACES arcs without RA/Dec values.  Ditto for GNIRS darks.
        ra = metadata_lookup.ra(storage_name.file_uri)
        dec = metadata_lookup.dec(storage_name.file_uri)
        if (ra is None and dec is None) or mode == 'imaging' or obs_type in ['ARC', 'BIAS', 'FLAT']:
            # DB 23-04-19
            # Ignore spatial WCS for the GRACES dataset with EPOCH=0.0.  Not important for a bias. For GMOS we skip
            # spatial WCS for biases  (and maybe for some other instruments).
            result = GracesSpectralTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
        else:
            result = GracesSpatialSpectralTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
    elif inst is Inst.MAROONX:
        obs_type = get_obs_type(metadata_lookup, storage_name)
        if obs_type == 'DARK':
            result = MAROONXTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
        elif obs_type in ['ARC', 'FLAT']:
            result = MAROONXSpectralTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
        else:
            result = MAROONXSpatialSpectralTemporal(
                storage_name, metadata_lookup, inst, clients, reporter, observation, config, provenance_finder
            )
    else:
        reporter.observable.rejected.record(mc.Rejected.MYSTERY_VALUE, storage_name.file_name)
        raise mc.CadcException(f'Mystery name {inst}.')
    logging.debug(f'Created {result.__class__.__name__} for mapping.')
    return [result]
