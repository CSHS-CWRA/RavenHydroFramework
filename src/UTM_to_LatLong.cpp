/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"

/* Ellipsoid model constants (actual values here are for WGS84) */
const double sm_a           = 6378137.0;         ///< Ellipsoid model constant for WGS 1984
const double sm_b           = 6356752.314;       ///< Ellipsoid model constant for WGS 1984
//const double sm_EccSquared  = 6.69437999013e-03; ///< Ellipsoid model constant for WGS 1984
const double UTMScaleFactor = 0.9996;            ///< UTM scale factor for WGS 1984

///////////////////////////////////////////////////////////////////
/// \brief Computes ellipsoidal distance from the equator to a point at a given latitude
/// \remark Translated from javascript of Chuck Taylor http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html \cite Taylor
/// \ref Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J., GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994. \cite Hoffmann-Wellenhof1994
/// \param phi [in] Latitude of point [rad]
/// \return Ellipsoidal distance of point from equator [m]
//
double ArcLengthOfMeridian (const double &phi)
{
  double alpha, beta, gamma, delta, epsilon, n;

  n = (sm_a - sm_b) / (sm_a + sm_b);
  alpha         = ((sm_a + sm_b) / 2.0)* (1.0 + (pow(n, 2.0) / 4.0) + (pow(n, 4.0) / 64.0));
  beta          = (-3.0 * n / 2.0) +
    (  9.0 * pow (n, 3.0) / 16.0 ) + ( -3.0 * pow (n, 5.0) / 32.0 );
  gamma         = ( 15.0 * pow (n, 2.0) / 16.0 ) + (-15.0 * pow (n, 4.0) / 32.0 );
  delta         = (-35.0 * pow (n, 3.0) / 48.0 ) + (105.0 * pow (n, 5.0) / 256.0);
  epsilon = (315.0 * pow (n, 4.0) / 512.0);

  return alpha
    * (phi + (beta * sin (2.0 * phi))
       + (gamma * sin (4.0 * phi))
       + (delta * sin (6.0 * phi))
       + (epsilon * sin (8.0 * phi)));
}

///////////////////////////////////////////////////////////////////
/// \brief Determines the central meridian for the given UTM zone
/// \remark Translated from javascript of Chuck Taylor http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
/// \remark Range of central meridian is the radian equivalent of [-177,+177]
/// \param zone [in] Integer value designating the UTM zone [1..60]
/// \return Central meridian for given UTM zone [rad], or zero if UTM zone is out of range
//
double UTMCentralMeridian (const int zone)
{
  return (-183.0 + (zone * 6.0))*DEGREES_TO_RADIANS;
}

///////////////////////////////////////////////////////////////////
/// \brief Computes the footpoint latitude for use in converting transver Mercator coordinates to ellipsoidal coordinates
/// \remark Translated from javascript of Chuck Taylor http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
/// \ref Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
/// \ref GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994. eqns. 10.18-10.22 \cite Hoffmann-Wellenhof1994
/// \param &y [in] The UTM northing latitude [m]
/// \return Footpoint latitude [rad]
//
double FootpointLatitude (const double &y)
{
  double y_, alpha_, beta_, gamma_, delta_, epsilon_, n;

  n = (sm_a - sm_b) / (sm_a + sm_b);// (Eq. 10.18)
  alpha_ = ((sm_a + sm_b) / 2.0) * (1 + (pow(n, 2.0) / 4) + (pow(n, 4.0) / 64));
  y_                    = y / alpha_;
  beta_         = (   3.0 *          n  / 2.0  ) + ( -27.0 * pow(n, 3.0) / 32.0) + (269.0 * pow(n, 5.0) / 512.0);
  gamma_        = (  21.0 * pow(n, 2.0) / 16.0 ) + ( -55.0 * pow(n, 4.0) / 32.0);
  delta_        = ( 151.0 * pow(n, 3.0) / 96.0 ) + (-417.0 * pow(n, 5.0) / 128.0);
  epsilon_= (1097.0 * pow(n, 4.0) / 512.0);

  return y_ + (beta_ * sin (2.0 * y_))
    + (gamma_ * sin (4.0 * y_))
    + (delta_ * sin (6.0 * y_))
    + (epsilon_ * sin (8.0 * y_));
}

///////////////////////////////////////////////////////////////////
/// \brief Converts a latitude/longitude pair to (x,y) coordintes in the Transverse Mercator (not UTM!) projection
/// \remark Translated from javascript of Chuck Taylor http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
/// \ref Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J., GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994. eqns. 10.18-10.22
///
/// \param phi [in] Latitude of point [rad]
/// \param lambda [in] Longitude of point [rad]
/// \param lambda0 [in] Central meridian [rad]
/// \param &x [out] Output x coordinate of point (transverse mercator)
/// \param &y [out] Output y coordinate of point (transverse mercator)
//
void MapLatLonToXY (const double phi,    //latitude, radians
                    const double lambda, //longitude, radians
                    const double lambda0,//central meridian, radians
                    double &x,
                    double &y)
{
  double N, nu2, ep2, t, t2, l;
  double l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;

  ep2 = (pow(sm_a, 2.0) - pow(sm_b, 2.0)) / pow(sm_b, 2.0);
  nu2 = ep2 * pow(cos (phi), 2.0);
  N = pow(sm_a, 2.0) / (sm_b * sqrt (1 + nu2));
  t = tan (phi);
  t2 = t * t;

  l = lambda - lambda0;

  /* Precalculate coefficients for l**n in the equations below
     so a normal human being can read the expressions for easting
     and northing
     -- l**1 and l**2 have coefficients of 1.0 */
  l3coef = 1.0 - t2 + nu2;
  l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
  l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2- 58.0 * t2 * nu2;
  l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2- 330.0 * t2 * nu2;
  l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
  l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);

  /* Calculate easting (x) */
  x = N * cos (phi) * l
    + (N / 6.0 * pow(cos (phi), 3.0) * l3coef * pow(l, 3.0))
    + (N / 120.0 * pow(cos (phi), 5.0) * l5coef * pow(l, 5.0))
    + (N / 5040.0 * pow(cos (phi), 7.0) * l7coef * pow(l, 7.0));

  /* Calculate northing (y) */
  y = ArcLengthOfMeridian (phi)
    + (t / 2.0 * N * pow(cos (phi), 2.0) * pow(l, 2.0))
    + (t / 24.0 * N * pow(cos (phi), 4.0) * l4coef * pow(l, 4.0))
    + (t / 720.0 * N * pow(cos (phi), 6.0) * l6coef * pow(l, 6.0))
    + (t / 40320.0 * N * pow(cos (phi), 8.0) * l8coef * pow(l, 8.0));
}

///////////////////////////////////////////////////////////////////
/// \brief Converts (x,y) coordintes in the Transverse Mercator (not UTM!) to projection a latitude/longitude pair
/// \remark Translated from javascript of Chuck Taylor http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
/// \ref Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
/// \ref GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994. eqns. 10.18-10.22 \cite Hoffmann-Wellenhof1994
///
/// \param &x [in] Output x ordinate of point (transverse mercator)
/// \param &y [in] Output y ordinate of point (transverse mercator)
/// \param lambda0 [in] Central meridian [rad]
/// \param phi [out] Latitude of point [rad]
/// \param lambda [out] Longitude of point [rad]
//
void MapXYToLatLon (const double &x,
                    const double &y,
                    const double &lambda0,
                    double &phi,
                    double &lambda)
{
  double phif, Nf, Nfpow , nuf2, ep2, tf, tf2, tf4, cf;
  double x1frac, x2frac, x3frac, x4frac, x5frac, x6frac, x7frac, x8frac;
  double x2poly, x3poly, x4poly, x5poly, x6poly, x7poly, x8poly;

  /* Get the value of phif, the footpoint latitude. */
  phif = FootpointLatitude (y);
  ep2 = (pow(sm_a, 2.0) - pow(sm_b, 2.0))/ pow(sm_b, 2.0);
  cf = cos (phif);
  nuf2 = ep2 * pow(cf, 2.0);
  Nf = pow(sm_a, 2.0) / (sm_b * sqrt (1 + nuf2));
  Nfpow  = Nf;
  tf = tan (phif);
  tf2 = tf * tf;
  tf4 = tf2 * tf2;

  /* Precalculate fractional coefficients for x^n in the equations
     below to simplify the expressions for latitude and longitude. */
  x1frac = 1.0 / (Nfpow  * cf);
  Nfpow  *= Nf;
  x2frac = tf / (2.0 * Nfpow );
  Nfpow  *= Nf;
  x3frac = 1.0 / (6.0 * Nfpow  * cf);
  Nfpow  *= Nf;
  x4frac = tf / (24.0 * Nfpow );
  Nfpow  *= Nf;
  x5frac = 1.0 / (120.0 * Nfpow  * cf);
  Nfpow  *= Nf;
  x6frac = tf / (720.0 * Nfpow );
  Nfpow  *= Nf;
  x7frac = 1.0 / (5040.0 * Nfpow  * cf);
  Nfpow  *= Nf;   // now equals Nf^8)
  x8frac = tf / (40320.0 * Nfpow );

  /* Precalculate polynomial coefficients for x**n.
     -- x**1 does not have a polynomial coefficient. */
  x2poly = -1.0 - nuf2;
  x3poly = -1.0 - 2 * tf2 - nuf2;
  x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2
    - 3.0 * (nuf2 *nuf2) - 9.0 * tf2 * (nuf2 * nuf2);
  x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;
  x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2+ 162.0 * tf2 * nuf2;
  x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);
  x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);

  /* Calculate latitude */
  phi = phif + x2frac * x2poly * (x * x)
    + x4frac * x4poly * pow(x, 4.0)
    + x6frac * x6poly * pow(x, 6.0)
    + x8frac * x8poly * pow(x, 8.0);

  /* Calculate longitude */
  lambda = lambda0 + x1frac * x
    + x3frac * x3poly * pow(x, 3.0)
    + x5frac * x5poly * pow(x, 5.0)
    + x7frac * x7poly * pow(x, 7.0);
}

///////////////////////////////////////////////////////////////////
/// \brief Converts a latitude/longitude pair to x,y coordinates in the UTM system
/// \remark Translated from Javascript of Chuck Taylor
/// http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
///
/// \param lat [in] Latitude, in decimal degrees
/// \param lon [in] Longitude, in decimal degrees
/// \param zone [in] UTM zone in which the point lies
/// \param &x [out] x-cooordinate in UTM system
/// \param &y [out] y-coordinate in UTM system
//
void LatLonToUTMXY (const double lat, //latitude, in decimal degrees
                    const double lon, //longitude, in decimal degrees
                    const int    zone,//UTM zone
                    double &x,
                    double &y)
{
  MapLatLonToXY (lat*DEGREES_TO_RADIANS, lon*DEGREES_TO_RADIANS, UTMCentralMeridian (zone), x,y);
  /* Adjust easting and northing for UTM system. */
  x = x * UTMScaleFactor + 500000.0;
  y = y * UTMScaleFactor;
  if (y < 0.0){ y = y + 10000000.0; }
}

///////////////////////////////////////////////////////////////////
/// \brief Converts (x,y) coordinates in UTM projection to a lat/long pair
/// \remark Translated from Javascript of Chuck Taylor
/// http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
///
/// \param &x [in] x-coordinate in UTM system
/// \param &y [in] y-coordinate in UTM system
/// \param zone [in] UTM zone in which the point lies
/// \param &southhemi [in] True if point is in southern hemisphere
/// \param lat [out] Latitude, in decimal degrees
/// \param lon [out] Longitude, in decimal degrees
//
void UTMXYToLatLon (const double &x,
                    const double &y,
                    const int    &zone,
                    const bool   &southhemi,
                    double &lat,
                    double &lon)
{
  double tmpx,tmpy;

  tmpx =x- 500000.0;
  tmpx /= UTMScaleFactor;

  /* If in southern hemisphere, adjust y accordingly. */
  tmpy=y;
  if (southhemi){tmpy -= 10000000.0;}
  tmpy /= UTMScaleFactor;

  MapXYToLatLon (tmpx, tmpy, UTMCentralMeridian(zone), lat,lon);
  lat*=RADIANS_TO_DEGREES;
  lon*=RADIANS_TO_DEGREES;
}
