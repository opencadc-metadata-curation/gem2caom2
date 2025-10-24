# How To Run gem2caom2

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. In the master branch of this repository, find the scripts directory, and copy the file gem_run_public.sh to the working directory. e.g.:

  ```
  wget https://raw.github.com/opencadc-metadata-curation/gem2caom2/master/scripts/gem_run_public.sh
  ```

2. Ensure the script is executable:

```
chmod +x gem_run_public.sh
```

3. To run the application:

```
./gem_run_public.sh
```
Running this will pull newly public FITS files and preview images from Gemini, and will generate thumbnail images from the preview.


# How To Run gem2caom2 Incrementally

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. In the master branch of this repository, find the scripts directory, and copy the file gem_run_incremental.sh to the working directory. e.g.:

  ```
  wget https://raw.github.com/opencadc-metadata-curation/gem2caom2/master/scripts/gem_run_incremental.sh
  ```

2. Ensure the script is executable:

```
chmod +x gem_run_incremental.sh
```

3. To run the application:

```
./gem_run_incremental.sh
```

Running this will generate CAOM observations for newly observed data at archive.gemini.edu.

# gem2caom2
An application to generate CAOM2 Observations from FITS file headers.

1. GMOS-N/S long-slit spectra:
   1. object: GN-2009B-Q-121-15-003
   1. arc:    GN-2010A-Q-35-10-002
   1. flat:   GN-2009B-Q-121-15-001
   1. dark:   GN-2007B-Q-112-14-018
1. GMOS-N/S multi-object spectra:
   1. object: GS-2003B-Q-5-29-003 (no OBSCLASS)
   1. object: GS-2005B-Q-27-4-001
   1. arc:    GS-CAL20060125-1-002
   1. flat:   GS-2005B-Q-27-33-001
   1. dark:   GS-2008A-Q-41-26-013
   1. mask:   GN2009AQ021-04
1. GMOS-N/S Integral-field spectra:
   1. object: GS-2005B-Q-22-17-006
   1. arc:    GS-2005B-Q-22-29-004
   1. flat:   GS-2005B-Q-54-2-004
   1. flat:   GS-2006A-Q-48-15-002 (twilight flat)
   1. dark:   GS-2009A-Q-30-6-007
1. NIRI imaging:
   1. object: GN-2009B-Q-52-332-011 (moving target; with AO so different FOV)
   1. object: GN-2011B-Q-28-10-002 (acquisition OBSCLASS)
1. NIRI long-slit spectra:
   1. object: GN-2008B-C-3-2-012
   1. flat:   GN-2008B-C-3-2-009
   1. arc:    GN-2010A-Q-44-175-015
1. GNIRS has no WCS info in extension; use primary header
1. GNIRS imaging:
   1. object: GN-2015B-SV-101-1061-005 (with AO)
   1. object: GN-CAL20151213-6-002 (without AO)
   1. object: GN-2010B-SV-142-10-007 (acquisition)
   1. object: GN-2010B-Q-2-44-003 (acqCal)
   1. dark:   GN-CAL20160202-3-039
   1. flat:   GN-2011A-Q-53-42-007
1. GNIRS spectroscopy:
   1. object: GN-2011B-Q-7-193-010 (cross-dispersed, short-slit)
   1. object: GN-2011B-Q-63-101-006 (long-slit)
   1. arc:    GN-2011B-Q-68-116-013 (cross-dispersed, short-slit)
   1. arc:    GN-2011B-Q-63-126-005 (long-slit)
   1. flat:   GN-2011B-Q-7-193-044 (cross-dispersed, short-slit)
   1. flat:   GN-2011B-Q-63-130-002 (long-slit)
1. GRACES has only telescope RA/Dec in the header, no WCS info.  Kludge perhaps a 5" square FOV consisting of a single pixel centered on RA/Dec values (it's a fiber-fed instrument)
1. GRACES spectroscopy: (No OBSCLASS keyword)
   1. object: GN-2015B-Q-1-12-1003 (issue: file name from jsonsummary != file name from fullheader)
   1. bias:   GN-CAL20150604-1000-1072
   1. flat:   GN-CAL20150604-1000-1081
   1. arc:    GN-CAL20150807-1000-1035
1. NIFS has no WCS info in extension; use primary header
1. NIFS integral field spectroscopy:
   1. object: GN-2012B-Q-51-108-001 (with AO)
   1. object: GN-2012B-Q-73-172-002 (without AO)
   1. dark:   GN-2016B-Q-27-31-007
   1. flat:   GN-2012A-Q-57-151-005
   1. flat:   GN-2011B-Q-59-70-016 (Ronchi flat)
   1. arc:    GN-2012A-Q-2-22-007
   1. image:  GN-2006A-C-11-670-006
1. GSAOI imaging:
   1. object: GS-2013B-Q-61-8-008
   1. object: GS-2012B-SV-499-21-002
   1. object: GS-2013B-DD-1-13-002 (acquisition)
   1. object: GS-CAL20181023-5-001 (subrasters; good test for spatial WCS)
   1. object: GS-CAL20140109-3-009 (domeflat; only identified by OBJECT keyword)
   1. flat:   GS-CAL20130201-3-017
   1. dark:   GS-2013B-Q-26-19-004
1. F2 (Flamingos 2) imaging:
   1. object: GS-2014B-Q-60-8-071
   1. object: GS-2015A-Q-60-126-001 (acquisition)
   1. object: GS-2014B-Q-17-69-002 (acqCal)
   1. dark:   GS-2014B-Q-17-53-030
   1. flat:   GS-2014B-Q-60-11-012
1. F2 long-slit spectroscopy:
   1. object: GS-2017A-Q-28-280-012
   1. object: GS-2016A-FT-18-54-002 (acqCal; very rare)
   1. object: GS-2016A-C-1-86-003 (acq; very rare)
   1. arc:    GS-2017B-Q-18-96-006
   1. flat:   GS-2017B-Q-45-156-079
1. F2 multi-object spectroscopy (rarely used apparently):
   1. object: GS-2015A-Q-63-328-014 (moving target; all seem to be!!)
   1. object: GS-2018B-SV-301-144-020
   1. flat:   GS-2015A-Q-63-365-016
   1. arc:    GS-2015A-Q-63-406-001
1. Note for GPI, below, datasets have two extensions. First is science image (with WCS), second is data quality (so datatype = Auxiliary) for each pixel (no WCS).
1. GPI polarimetric imaging:
   1. object: GS-2018A-FT-101-5-043
   1. dark:   GS-CAL20150410-1-007
   1. flat:   GS-CAL20140323-6-014
   1. arc:    GS-CAL20160224-12-017
1. GPI integral field spectra:
   1. object: GS-2014A-SV-408-6-003
   1. dark:   GS-CAL20140320-7-011
   1. flat:   GS-CAL20140315-4-004
   1. arc:    GS-CAL20140316-5-012
1. The challenge for NICI? Two extensions contain data but with different filters!
1. NICI (retired) imaging:
   1. object: GS-2009B-Q-14-129-029
   1. flat:   GS-CAL20100110-1-026
   1. dark:   GS-2013A-Q-39-158-008
1. Michelle is a retired visitor instrument.  Spatial WCS info is in primary header.  There are a variable number of FITS extensions defined by primary keyword NUMEXT; assume the same spatial WCS for each for now - it differs only slightly because of telescope 'chopping' and 'nodding' used in acquisition. 
1. Question:  Do we 'edit' the instrument name?  e.g. instead of 'michelle' use 'Michelle'?
I guess not since Gemini displays result set with 'michelle' from the header.
1. Michelle imaging:
   1. object: GN-2008A-Q-43-6-002
   1. object: GN-2009B-C-1-62-001 (acquisition)
   1. flat:   GN-CAL20060413-7-001 (rare)
   1. bias:   GN-CAL20080308-8-016
1. Michelle long-slit spectroscopy:
   1. object: GN-2006A-Q-58-9-007
   1. object: GN-2007A-C-11-234-001 (acquisition)
   1. flat:   GN-2010B-C-3-21-002
   1. bias:   GN-2006A-C-14-49-002
1. T-ReCS is that, like Michelle, it has a variable number of extensions depending on the amount of chopping/nodding cycles carried out during the observation.  Use primary header for this one for the WCS info.  It can be assumed the same for each extension.  The maximum chop throw for Gemini is 15 arcsec I believe so not worth trying to adjust each extensions WCS solution by this amount (it likely can’t even be done with info in the headers).
1. T-ReCS imaging:
   1. object: GS-CAL20050102-1-001 (before OBSCLASS keyword)
   1. object: GS-2005A-Q-15-1-001
   1. object: GS-2005B-Q-10-22-003 (acquisition)
1. T-ReCS spectroscopy:
   1. object: GS-2012A-Q-7-31-001
   1. object: GS-2008A-C-5-35-002 (acquisition)
1. bHROS spectroscopy:
   1. object: GS-2006B-Q-47-76-003
   1. flat:   GS-2006B-Q-7-32-008
   1. arc:    GS-2005B-SV-302-20-001
   1. bias:   GS-2005B-SV-301-16-005
1. Question:  Do we 'edit' the instrument name?  e.g. instead of 'hrwfs' use 'HRWFS'?
Again, likely not since Gemini displays results with 'hrwfs' from the header.
1. hrwfs/Acquisition Camera imaging:
   1. object: GS-2003A-Q-6-1-002
   1. bias:   GS-CAL20030303-7-0006 (old style filename/header)
   1. bias:   GS-CAL20030730-10-006 (new style filename/header)
   1. flat:   GS-CAL20031218-1-034
   1. dark:   GS-CAL20030105-2-0072
1. Note: PHOENIX headers use non-standard OBS_TYPE keyword (instead of OBSTYPE).  And to make it worse, this is always set to 'Object'.  Sigh.  We could look for 'flat', 'dark', 'arc/comp/lamp' or 'comparison' in the OBJECT keyword, but this would be subject to error (at least if we assumed all others were 'object' exposures.
1. PHOENIX spectroscopy:
   1. object: GS-2006B-C-8-2-0052
   1. flat:   GS-2006A-DD-1-1-0073
   1. arc:    GS-2003B-Q-51-27-0073
   1. dark:   GS-2006A-C-10-1-0258
   1. GS-2002A-DD-1-17-0171:  DARK because VIEW_POS contains ‘dark’, even though ‘object’ is ‘gcal 6420’
   1. GS-CAL20020512-9-0277:  also a dark because json ‘object’ is ‘dark (gcal4618)’.
1. Very early instruments:
   1. Hokupa'a + QUIRC, OSCIR and Flamignos were among the first instruments used on Gemini.  All were visitor instruments.  For OSCIR and Hokupa'a + QUIRC there is no spatial WCS in the headers.  For OSCIR the (non-standard) OBS_TYPE keyword is always 'object' and the filter name appears to be frequently missing. Flamingos also suffers from the incorrect use of the non-standard OBS_TYPE keyword, and Hokupa'a + QUIRC often does as well. For the latter they use IMAGETYP as the keyword but I see only 'dark' and 'object' values.  Gemini metadata also does not tell us what the observing mode is for any of these (imaging vs. spectroscopy).  
   1. Note: it is always imaging for Hokupa'a + QUIRC so that could be hardcoded.  I've given only one sample file for OSCIR/Flamingos since we can't distinguish between different observation types.  Gemini metadata does NOT include the filter name for any of these instruments, although this appears to be available in the Hokupa'a+QUIRC 'FILTER' keyword.  
   1. Since FILTER is in the Hokupa'a + QUIRC header you might be able to use the SVO's 'generic' filters to determine energy information:  http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=Generic&gname2=Johnson_UBVRIJHKL
1. OSCIR imaging/spectroscopy??:
   1. object: GS-2001B-Q-31-9-007
1. Flamingos Imaging
   Flamingos Imaging:
   1. object: GS-2002B-DD-2-6-0230
   1. object: GS-2002B-Q-5-20-0205
   1. object: GS-2002A-Q-13-2-0186
   1. flat:     GS-CAL20020707-8-0034
   1. dark:    GS-CAL20020707-9-0040
   1. object: GS-2002B-DD-2-6-0230
1. FLAMINGOS long-slit spectra:
   1. object: GS-2002A-Q-7-1-0071
   1. object: GS-2002B-DS-1-46-0031 (sample acquisition observation -> image)
   1. object: GS-2002B-DS-1-46-0042 (resulting spectrum)
   1. arc:    GS-CAL20020625-11-0121
   1. flat:    GS-CAL20020625-10-0115
1. FLAMINGOS multi-object spectra:
   1. object: GS-2002B-DS-1-31-0072
   1. object: GS-2002B-DS-1-41-0011 (sample acquisition observation -> image)
   1. object: GS-2002B-DS-1-41-0016 (spectrum of above standard star)
   1. flat: GS-CAL20021105-1-0001
   1. flat: GS-CAL20021112-3-0024
1. Hokupa'a + QUIRC imaging:
   1. object: GN-2002A-DD-1-319-591
   1. dark:   GN-CAL20020424-1-007
1. TEXES was a more recent visitor instrument so should be more consistent.
1. TEXES
   1. flat field: TX20131117_flt.3002.fits
   1. unprocessed science: TX20131117_raw.3002.fits
   1. flat field: TX20170321_flt.2507.fits
   1. unprocessed science: TX20170321_raw.2507.fits
   1. reduced science: TX20170321_red.2507.fits
1. ABU returns no data on the Gemini site.
1. ABU - TBD
1. A quick look suggest CIRPASS will be the same as the other older visitor instruments.   e.g. no use of IMAGTYPE or OBSTYPE keyword, no WCS spatial information at all as for OSCIR.  
1. CIRPASS
   1. GS-CAL20030320-3-2731    dome flat
   1. GS-CAL20030630-4-3463    another dome flat
   1. GS-CAL20030309-9-1247     dark
   1. GS-2003A-Q-10-19-1204      science object
   1. GS-2003A-Q-3-22-3598         science object
   1. GS-CAL20030315-18-2027   standard star (obstype=OBJECT but intent still CALIBRATION)
   1. GS-CAL20030630-1-3384      flat
   1. GS-CAL20030313-3-1769      arc
1. GMOS Processed Files, Separate Composite Observations
   1. GS-CAL20161227-5-001                S20161227S0051
   1. GS-CAL20161227-5-001-RG-FRINGE      rgS20161227S0051_fringe
   1. GS-CAL20131007-900-067-FRINGE       S20131007S0067_fringe
   1. GN-CAL20120320-900-328-STACK-FRINGE N20120320S0328_stack_fringe
   1. GN-CAL20110927-900-170-FRINGE       N20110927S0170_fringe
   1. GS-CAL20181016-5-001-MRG-FRINGE     mrgS20181016S0184_fringe
   1. GS-2016B-Q-72-23-001-MRG-ADD        mrgS20160901S0122_add
   1. GN-2016A-Q-68-46-001-MRG-ADD        mrgN20160311S0691_add
   1. GS-2016A-Q-7-175-001-MFRG-ADD       mfrgS20160310S0154_add
   1. GN-2005B-Q-28-32-001-MRG-ADD        mrgN20050831S0770_add
   1. GS-2004B-Q-42-1-001-MFRG-ADD        mfrgS20041117S0073_add
   1. GN-2002A-SV-78-7-003-FMRG-ADD       fmrgN20020413S0120_add
   1. GS-2012B-Q-1-32-002-MRG             mrgS20120922S0406
   1. GS-2012B-Q-1-32-002                 S20120922S0406
   1. GN-2004B-Q-30-1-001-MRG             mrgN20041016S0095
   1. GN-CAL20160404-7-017-FLAT-PASTED    N20160403S0236_flat_pasted
   1. GS-CAL20131109-17-001-RG-FRINGE     rgS20131109S0166_FRINGE
   1. GS-CAL20141226-7-026-G-BIAS         GS20141226S0203_BIAS
   1. GN-CAL20141109-2-001-BIAS           N20141109S0266_bias
1. GMOS Processed Files, Additional Planes on existing Simple Observations
   1. GS-2010A-Q-36-6-358-RG              rgS20100316S0366
1. TReCS Processed Files, Additional Planes on existing Simple Observations
   1. GS-2012B-Q-90-366-003-R             rS20121030S0136 
1. NIRI Processed
   1. GN-2013B-Q-75-163-011_STACK         N20140313S0072_flat
   1. GN-2015B-Q-53-138-061_STACK         N20150804S0348_dark
   1. GN-2007B-Q-107-150-004_DARK         N20070819S0339_dark
   1. GN-2013A-Q-63-54-051_FLAT           N20130404S0512_flat
1. Flamingos Processed
   1. GS-2013B-Q-16-277-019_STACK         S20140124S0039_dark
   1. GS-CAL20141129-1-001_DARK           S20141129S0331_dark
   
