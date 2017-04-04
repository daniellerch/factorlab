/*
   Copyright 2015 Daniel Lerch Hostalot.

   This file is part of factorlab.

   factorlab is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   factorlab is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with factorlab; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */



#ifndef __ECM_HPP__
#define __ECM_HPP__

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>


NTL::ZZ ecm(NTL::ZZ &n, long B1, long B2);
NTL::ZZ ecrho(NTL::ZZ &n, long B1, long B2);


#endif


