#pragma once

#include <array>
#include <limits>
#include <stdexcept>

/*
	http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
    https://oeis.org/A039598
	https://oeis.org/A039598/b039598.txt

    D = [0]*(n+2); D[1] = 1
    b = True; h = 1
    for i in range(2*n) :
        if b :
            for k in range(h, 0, -1) : D[k] += D[k-1]
            h += 1
        else :
            for k in range(1, h, 1) : D[k] += D[k+1]
        b = not b
        if b : print([D[z] for z in (1..h-1) ])
*/

/*
351 69533550916004
352 124680849918352
353 155851062397940
354 160878516023680
355 144539291740025
356 115631433392020
357 83322650532485
358 54414792184480
359 32308782859535
360 17464206951100
361 8594228157515
362 3846367846720
363 1562586937730
364 574609830760
365 190559382650
366 56724653440
367 15067486070
368 3545290840
369 732179630
370 131185600
371 20087795
372 2576860
373 269399
374 22048
375 1325
376 52
377 1
 */

template<typename T, unsigned n>
constexpr std::array<T,n> calc_last_line_of_catalan_triangle()
{
  std::array<T,n+2> D{0};

  bool b = true;
  int h = 1;
  D.at(1) = 1;

  for( int i = 0; i < 2 * int(n); ++i ) {
    if( b ) {
      for( int k = h; k > 0; --k ) {

        if( std::numeric_limits<T>::max() - D.at(k-1) < D.at(k) ) {
          throw std::overflow_error("too much");
        }

        D.at(k) += D.at(k-1);
      }
      h += 1;
    } else {
      for( int k = 1; k < h; ++k ) {

        if( std::numeric_limits<T>::max() - D.at(k+1) < D.at(k) ) {
          throw std::overflow_error("too much");
        }

        D.at(k) += D.at(k+1);
      }
    }

    b = !b;
  }

  std::array<T,n> ret{0};

  for( int z = 1,i = 0; z < h; ++z, ++i ) {
    ret.at(i) = D.at(z);
  }

  return ret;
}
