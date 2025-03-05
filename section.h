/*! \file section.h
    \brief Definition of Poincare sections.

    $Author: roldan $
    $Date: 2013-03-11 11:36:04 $
*/

#ifndef SECTION_H_INCLUDED
#define SECTION_H_INCLUDED

/** A Poincare section for the RTBP is given by a plane (with orientation), and
 * consists in a base point p(x,y,z) and a unit normal vector n(x,y,z) to the
 * plane. 
 */
typedef struct {
	double p[3];	/**< base point */
	double n[3];	/**< normal vector */
} section_t;

#endif // SECTION_H_INCLUDED
