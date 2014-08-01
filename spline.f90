## Copyright (C) 1996, 1997 John W. Eaton
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.

## -*- texinfo -*-
##  at deftypefn {Function File} {@var{pp} =} spline (@var{x}, @var{y})
##  at deftypefnx {Function File} {@var{yi} =} spline (@var{x}, @var{y}, @var{xi})
## Piece-wise interpolation with cubic splines.
##
## The argument  at var{x} should be a vector and @var{y} should either be
## a vector of the same length as  at var{x} or a matrix with length
## ( at var{x}) number of columns. The spline-function uses cubic spline
## interpolation to find  at var{yi} corresponding to @var{xi}. If the
## third argument is omitted, spline returns a piece-wise polynom on
## pp-form.
##
## Example:
## 
##  at example
## x  = -1:0.25:1;
## xi = -1:0.01:1;
## y  = [exp(x)./(1 + 9*x.^2); sin(4*x)];
## y1 = [exp(xi)./(1 + 9*xi.^2); sin(4*xi)];
## yi = spline (x, y, xi);
## plot (xi, y1, ["g-;calculated;";"g-;;"], \
##       xi, yi, ["b-;interpolated;";"b-;;"], \
##       x, y, ["ro;interpolation-points;";"ro;;"])
##  at end example
##  at end deftypefn
##  at seealso{ppval}

## Author: Mats Jansson <jansson at ieee dot org>
## Created: March 5, 2000 - April 6, 2001
## Keywords: piecewise interpolation cubic spline

function yi = spline (x, y, xi)

  if (nargin < 2),
    usage ("spline: expecting 2 or 3 arguments.");
  endif
  
  x = x(:); ## Make sure x is a column vector
  n = size (x, 1);
  if (n < 2),
    usage ("yi = spline (x, y, xi): x should have at least 2 elements.");
  endif
  
  if (is_vector (y)),
    y = y(:);    ## Make sure y is a columnvector
    [nry, ncy] = size (y);
    if (nry != n),
      error ("yi = spline (x, y, xi): x and y should have the same number of elements.");
    endif
  elseif (is_matrix (y)),
    [nry, ncy] = size (y);
    if (ncy != n),
      ## The number of columns of y should match the length of x (I
      ## think this is the way Matlabs' spline-function behaves).
      error ("yi = spline (x, y, xi): The number of columns of y should match the length of x.");
    else
      y = y.';
      [nry, ncy] = size (y);
    endif
  endif
  
  
  [x, ind] = sort (x);
  if (! all (diff (x))),
    usage ("yi = spline (x, y, xi): x should be distinct.");
  endif
  
  y = y(ind,:);
  
  h = diff (x);
  d = diff (y);
  
  if (n == 2),
    ## Linear interpolation
    coef = [d./h(:,ones (1, ncy)); y(1,:)].';
  else
    ## We want the differentiate and the second differentiate to be
    ## continuous in each point. To find the differentiate, k, in each
    ## point, we have to solve the equation system:
    ##
    ## H * k = 3b, where H is the n x n tri-diagonal matrix
    ##   _                                                           _
    ##  |  2h(1)       h(1)          0         0      ...        0    |
    ##  |   h(2)   2h((2)+h(1))     h(1)       0      ...        0    |
    ##  |    0         h(3)     2(h(3)+h(2))  h(2)    ...        0    |
    ##  |    .          .               .                        .    |
    ##  |    .          .                  .                     .    |
    ##  |    .          .         h(n-1)  2(h(n-1)+h(n-2)    h(n-2)   |
    ##  |_   0          0           0          h(n-1)       2h(n-1)  _|
    ##
    ## and b is the n-by-ncy matrix
    ##   _                                              _
    ##  |                     d(1,:)                     |
    ##  |        h(1)/h(2)*d(2,:)+h(2)/h(1)*d(1,:)       |
    ##  |        h(2)/h(3)*d(3,:)+h(3)/h(2)*d(2,:)       |
    ##  |                      ...                       |
    ##  |  h(n-2)/h(n-1)*d(n-1,:)+h(n-1)/h(n-2)*d(n-2,:) |
    ##  |_                   d(n-1,:)                   _|
    ##
    
    H = diag (2 * [h(1); h(2:n-1)+h(1:n-2); h(n-1)]) + \
        diag ([h(1); h(1:n-2)], 1) + \
        diag ([h(2:n-1); h(n-1)], -1);

    b = zeros (n, ncy);
    
    ## Repeat the column-vector h ncy times to form a matrix of the
    ## same size as d.
    h = h(:,ones (1, ncy));

    b(1,:) = d(1,:);
    b(n,:) = d(n-1,:);
    b(2:n-1,:) = h(1:n-2,:)./h(2:n-1,:).*d(2:n-1,:) + \
        h(2:n-1,:)./h(1:n-2,:).*d(1:n-2,:);
    
    ## H is a sparse matrix, this system could be solved more efficient.
    k = H \ (3*b);
    
    ## We should have n-1 3:rd degree interpolation polynoms for each
    ## column of y.
    coef = \
        [((h(1:n-1,:).*(k(1:n-1,:)+k(2:n,:))-2*d(1:n-1,:)) ./ \
          h(1:n-1,:).^3).'(:); \
         ((3*d(1:n-1,:)-2*h(1:n-1,:).*k(1:n-1,:)-h(1:n-1,:).*k(2:n,:)) ./ \
          h(1:n-1,:).^2).'(:); \
         k(1:n-1,:).'(:); \
         y(1:n-1,:).'(:)];
  endif
  
  yi = mkpp (x, coef, ncy);

  if (nargin >= 3),
    yi = ppval (yi, xi);
  endif
endfunction


