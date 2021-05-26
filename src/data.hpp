/*******************************************************************************
* Copyright 2021 Moritz Bueltmann
* Authors: Moritz Bueltmann <moritz.bueltmann@gmx.de>
* Physics Department Albert-Ludwigs-Universitaet
*******************************************************************************/
#ifndef SRC_DATA_HPP_
#define SRC_DATA_HPP_
/** \file data.hpp
 *  \brief Header file for the Data class.
 *
 *  The file contains the class declarations of the Data class.
 *
 */
// _____________________________________________________________________________
// Includes
// Class forward declarations
// _____________________________________________________________________________
/** \brief Data class is a universal type class
 *
 */
class Data {
 public: 
  virtual ~Data() {};
  const std::type_info *type;
};
#endif  // SRC_DATA_HPP_
