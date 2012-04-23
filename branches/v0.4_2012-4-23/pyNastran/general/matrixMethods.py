## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
def gauss(nGauss):
    """
    http://en.wikipedia.org/wiki/Gaussian_quadrature
    """
    if nGauss=1:
        P = [0.]
        W = [2]
    elif nGuass==2:
    
        p = 1./sqrt(3)
        P = [w,-w]
    elif nGuass==3:
        p = sqrt(3/5.)
        P = [p,0.,-p]
    elif nGuass==4:
        p1 = (3-2.*sqrt(6/5))/7.
        p2 = (3+2.*sqrt(6/5))/7.
        P = [-p1,p1,-p2,p2]
        w1 = (18+sqrt(30))/36.
        w2 = (18-sqrt(30))/36.
        W = [w1,w1,w2,w2]
    elif nGuass==5:
        p1 = 1/3.*sqrt(5-2*sqrt(10./7.))
        p2 = 1/3.*sqrt(5+2*sqrt(10./7.))
        w1 = (322+13*sqrt(70))/900.
        w2 = (322-13*sqrt(70))/900.
        P = [0.,-p1,p1,-p2,p2]
        W = [128./225.,w1,w1,w2,w2]
    else:
        raise NotImplementedError('the code only supports up to 5 quadrature points')
    return (P,W)
