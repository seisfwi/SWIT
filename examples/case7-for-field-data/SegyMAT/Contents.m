% SegyMAT : A toolbox to read, write and manipulating SEG Y formatted files
% Version 1.00
%
% New Features.
%   README      
%
% Main
%   ReadSegy          - Reads Segy File 
%   ReadSegyHeader    - Reads SegyHeader from Segy File 
%   ReadSegyFast      - Reads Segy File in fast mode. No header values will be read.
%   WriteSegy         - Write Segy formatted data   
%   WriteSegyStructure- Write Segy formatted data using SegyMAT data structures
%
%   ReadSu            - Reads a SU formatted file. 
%   ReadSuFast        - Reads a SU formatted file in fast mode. No header values will be read.
%   WriteSu           - Write SU formatted data   
%   WriteSuStructure  - Write Su formatted data using SegyMAT data structures
%
% Lower Level IO
%   GetSegyHeader          - Reads the segyheader of a SEGY Y formatted file 
%   GetSegyHeaderBasics    - Default Segy header settings
%
%   GetSegyTrace.m         - Read Segy Trace Header and Data from filehandle
%   GetSegyTraceHeader     - Read Segy Trace Header from filehandle
%   GetSegyTraceData       - Read Segy Trace Data from filehandle
%
%   PutSegyHeader          - Write Segy Header to filehandle
%   PutSegyTrace           - Write Segy Trace Header and Data to filehandle
%
%   InitSegyTraceHeader    - Initalize all fields in the SegyTraceheader
%   CheckSegyTraceHeader.m - Check a SegyTraceHeader for all required fields
%
%
% SU<-> SEG-Y conversion
%   Su2Segy - Convert SU formatted files to SEG Y
%   Segy2Su - Convert SEG Y formatted files to SU
%
% Plotting
%   wiggle  - wiggle/variable area/image plotting of seismic data
%
% Misc
%   ibm2num - Convert IBM 32 bit floatto double
%   num2ibm - Convert IEEE 754 doubles to IBM 32 bit floating point format
%   ebcdic2ascii - convert ebcdic to ascii format
%   SegymatVerbose - controls amount of info written to screen
%   SegymatVersion - Return the current SegyMAT version
%
%
% Seismic Processing :
%   SegyMAT_GAIN : 'agc' and 'power' gain.
%
% (C) 2001-2004 Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
%

% 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%

% Demonstrations.
%   SegyMATdemo1 - Writes a Segy formatted filr, reads it and displays it
%   (NOT IMPLEMENTED)

