// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_PROPERTIES_HPP_
#define SRC_PROPERTIES_HPP_
/** \file properties.hpp
 *  \brief Header file for the Properties class.
 *
 *  The file contains the class declarations of the Properties class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <stdexcept>
#include <iostream>
#include <string>
#include <unordered_map>
#include <typeinfo>
// Class forward declarations
// _____________________________________________________________________________
/** \brief Properties class is a container for species or system properties
 *
 */
class Properties {
  /** \brief The System class has full access on the properties. 
   *
   *  In general, the properties must not been changed once they are defined to 
   *  avoid manipulation of calculations during runtime. However, the System 
   *  needs to have access on properties in order to, for instance, update the 
   *  valancies of all species. For this reason, System is declared as a friend.
   *
   */
  // TODO(Andreas): Keine gute idee, wenn System ein template ist.
  //                Warum nicht die update-Funktion public machen?
  // friend class System;
 public:
  /** \brief Constructors
   *
   */
  Properties();
  Properties(const Properties& other);
  /** \brief Destructor
   *
   */
  ~Properties();
  /** \brief Removes all properties
   *
   */
  void clear();
  /** \brief Returns true if property is contained
   *
   */
  bool contains_property(const std::string& property_name) const;
  /** \brief Add a property with an arbitrary data type
   *
   */
  template<typename T>
  void add_property(
      const std::string& property_name, T property_value) {
    TemplateData<T>* new_property = new TemplateData<T>(property_value);
    properties[property_name] = new_property;
  }
  /** \brief Returns a property with an arbitrary data type
   *
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value) const {
    if (contains_property(property_name)) {
      if (typeid(T) != *(properties.at(property_name)->type)) {
        std::cerr << "Properties::get_property(): \"";
        std::cerr << "Error: Requested type and property type not the same.\"";
        std::cerr << std::endl;
        exit(1);
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          properties.at(property_name)))->value;
      return true;
    } else {
      throw &missing_property_error;
      return false;
    }
  }
  /** \brief Data class is a universal type class. 
   *
   *  The Data class is used to derive a TemplateData class that stores data of
   *  a certain type (defined by the user via the template declaration). 
   *  For each TemplateData class is also a Data class, different realizations 
   *  (of different data type) of TemplateData can be stored in a container that
   *  holds Data objects. 
   *
   */
  class Data {
   public:
    virtual ~Data() {}
    const std::type_info *type;
  };
  /** \brief TemplateData class is a container to hold an object of arbitrary 
   *         type. 
   *
   *  TemplateData is a container to hold an object of arbitrary type. In order
   *  to allow to store several TemplateData realizations of different type in 
   *  one container (like a vector or map), TemplateData is a child of the Data
   *  class. 
   *
   */
  template <typename T> class TemplateData : public Data {
   public:
    /** \brief Constructor
     *
     */
    explicit TemplateData(T value) {
      type = &(typeid(T));
      this->value = value;
    }
    /** \brief Stores the value of template type
     *
     */
    T value;
  };
  /** \brief Contains all properties.
   *
   */
  std::unordered_map<std::string, Data*> properties;
  /** \brief std::exception MissingPropertyException.
   */
  class MissingPropertyException : public std::exception {
   public:
    /** \brief Overwrite the exception information function what().
     *
     *  \return the text "The adressed property does not exist.".
     */
    virtual const char* what(void) const throw() {
      return "The requested property does not exist.";
    }
  }
  /** \brief Exception MissingPropertyException missing_property_error.
   *
   *  The exception is thrown if a requested property is missing. 
   */
  missing_property_error;

 protected:
  /** \brief Update a property with an arbitrary data type. 
   * 
   *  If the property already exists, its value is updated. 
   *  Otherwise, the property is added.
   * 
   *  \param property_name Name of the property. 
   *  \param property_value Value of the property. 
   *
   *  Example for calling the function: <br>
   *  update_property<double>("my property", 4.5);
   *
   */
  template<typename T>
  void update_property(
      const std::string& property_name, T property_value) {
    auto search = properties.find(property_name);
    if (search != properties.end()) {
      // Property exists: update value
      delete(properties[property_name]);
      TemplateData<T>* new_property = new TemplateData<T>(property_value);
      properties[property_name] = new_property;
    } else {
      // Property does not exist: add value
      this->add_property<T>(property_name, property_value);
    }
  }
};
#endif  // SRC_PROPERTIES_HPP_
