% ascii2ebcdic : Converts ASCII formatted text to EBCDIC formatted text
%
% CALL : ebcdic=ascii2ebcdic(ascii);
%
% ascii  : Array on unsigned integers
% ebcdic : Array on unsigned integers
%
% (C) 2002-2009, Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas.mejer.hansen@gmail.com
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
function ebcdic=ascii2ebcdic(ascii)
  
  ebcdic_arr=[
    '00';'01';'02';'03';'04';'05';'06';'07';'08';'09';
    '0A';'0B';'0C';'0D';'0E';'0F';'10';'11';'12';'13';
    '14';'15';'16';'17';'18';'19';'1A';'1B';'1C';'1D';
    '1E';'1F';'20';'21';'22';'23';'24';'25';'26';'27';
    '28';'29';'2A';'2B';'2C';'2D';'2E';'2F';'30';'31';
    '32';'33';'34';'35';'36';'37';'38';'39';'3A';'3B';
    '3C';'3D';'3E';'3F';'40';'4A';'4B';'4C';'4D';'4E';
    '4F';'50';'5A';'5B';'5C';'5D';'5E';'5F';'60';'61';
    '6A';'6B';'6C';'6D';'6E';'6F';'79';'7A';'7B';'7C';
    '7D';'7E';'7F';'81';'82';'83';'84';'85';'86';'87';
    '88';'89';'91';'92';'93';'94';'95';'96';'97';'98';
    '99';'A1';'A2';'A3';'A4';'A5';'A6';'A7';'A8';'A9';
    'C0';'C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8';'C9';
    'D0';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';'D9';
    'E0';'E2';'E3';'E4';'E5';'E6';'E7';'E8';'E9';'F0';
    'F1';'F2';'F3';'F4';'F5';'F6';'F7';'F8';'F9';'FF';
    '00';'01';'02';'03';'37';'2D';'2E';'2F';'2F';'16';
    '05';'25';'0B';'0C';'0D';'10';'11';'12';'13';'3C';
    '3D';'32';'26';'18';'3F';'27';'1C';'1D';'1D';'1E';
    '1F';'07';'40';'5A';'7F';'7B';'5B';'6C';'50';'7D';
    '4D';'5D';'5C';'4E';'6B';'60';'4B';'61';'F0';'F1';
    'F2';'F3';'F4';'F5';'F6';'F7';'F8';'F9';'7A';'5E';
    '4C';'7E';'6E';'6F';'7C';'00';'E0';'00';'00';'6D';
    '79';'C0';'4F';'D0';'A1'];

ascii_arr=[
    '00';'01';'02';'03';'9C';'09';'86';'7F';'97';'8D';
    '8E';'0B';'0C';'0D';'0E';'0F';'10';'11';'12';'13';
    '9D';'85';'08';'87';'18';'19';'92';'8F';'1C';'1D';
    '1E';'1F';'80';'81';'82';'83';'84';'0A';'17';'1B';
    '88';'89';'8A';'8B';'8C';'05';'06';'07';'90';'91';
    '16';'93';'94';'95';'96';'04';'98';'99';'9A';'9B';
    '14';'15';'9E';'1A';'20';'A2';'2E';'3C';'28';'2B';
    '7C';'26';'21';'24';'2A';'29';'3B';'AC';'2D';'2F';
    'A6';'2C';'25';'5F';'3E';'3F';'60';'3A';'23';'40';
    '27';'3D';'22';'61';'62';'63';'64';'65';'66';'67';
    '68';'69';'6A';'6B';'6C';'6D';'6E';'6F';'70';'71';
    '72';'7E';'73';'74';'75';'76';'77';'78';'79';'7A';
    '7B';'41';'42';'43';'44';'45';'46';'47';'48';'49';
    '7D';'4A';'4B';'4C';'4D';'4E';'4F';'50';'51';'52';
    '5C';'53';'54';'55';'56';'57';'58';'59';'5A';'30';
    '31';'32';'33';'34';'35';'36';'37';'38';'39';'9F';
    '00';'01';'02';'03';'04';'05';'06';'07';'07';'08';
    '09';'0A';'0B';'0C';'0D';'10';'11';'12';'13';'14';
    '15';'16';'17';'18';'1A';'1B';'1C';'1D';'1D';'1E';
    '1F';'7F';'20';'21';'22';'23';'24';'25';'26';'27';
    '28';'29';'2A';'2B';'2C';'2D';'2E';'2F';'30';'31';
    '32';'33';'34';'35';'36';'37';'38';'39';'3A';'3B';
    '3C';'3D';'3E';'3F';'40';'5B';'5C';'5D';'5E';'5F';
    '60';'7B';'7C';'7D';'7E'
];

ascii_dec_array=hex2dec(ascii_arr);

ebcdic_dec_array=hex2dec(ebcdic_arr);


for i=1:length(ascii)
  m=find(ascii_dec_array==ascii(i));
  if length(m)>1, 
    m=m(1);
  end
  if length(m)==0,
    ebcdic(i)=0;
  else
    ebcdic(i)=(ebcdic_dec_array(m));  
  end
end
