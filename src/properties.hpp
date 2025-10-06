// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
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
#include "../../parameter_handler/src/parameter_handler.hpp"
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
  Properties(ParameterHandler* parameters);
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
  /** \brief Add an indexed property
   */
  template<typename T>
  void add_property(
      const std::string& property_name, T property_value, size_t index) {
    add_property(indexed(property_name, index), property_value);
  }
  /** \brief Returns a property with an arbitrary data type
   *
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value) const {
    if (contains_property(property_name)) {
      if (typeid(T) != *(properties.at(property_name)->type)) {
        std::cerr << "Properties::get_property(): \"";
        std::cerr << "Error: Requested type and property type not the same for ";
	std::cerr << "property '" << property_name << "'.\"";
        std::cerr << std::endl;
        exit(1);
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          properties.at(property_name)))->value;
      return true;
    } else {
      std::cout << "Properties: " << property_name << " not specified.\n";
      throw &missing_property_error;
      return false;
    }
  }
  /** \brief Returns a property, possibly extracting it from parameter_handler
   *
   * Note that the ParameterHandler usually needs parameter names without spaces
   * therefore this function also attempts to find the parameter by an
   * alternative name without spaces.
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value) {
    if (contains_property(property_name)) {
      // extract property from properties
      if (typeid(T) != *(properties.at(property_name)->type)) {
        std::cerr << "Properties::get_property(): \"";
        std::cerr << "Error: Requested type and property type not the same for ";
	std::cerr << "property '" << property_name << "'.\"";
        std::cerr << std::endl;
        exit(1);
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          properties.at(property_name)))->value;
      return true;
    } else if (parameter_handler != nullptr) {
      // extract property from parameter handler & add it to properties for
      // future reuse
      if (extract_property(property_name, property_value)) {
	add_property(property_name, *property_value);
        return true;
      // try the same with alternative name
      } else if (extract_property(parameter_name(property_name),
			          property_value)) {
	add_property(property_name, *property_value);
        return true;
      }
    }
    // else:
    std::cout << "Properties: " << property_name << " not specified.\n";
    throw &missing_property_error;
    return false;
  }
  /** \brief Returns an indexed property
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value,
		  size_t index) const {
    return get_property(indexed(property_name, index), property_value);
  }
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value,
		  size_t index) {
    return get_property(indexed(property_name, index), property_value);
  }
  /** \brief Specifies labeling convention for indexed properties
   */
  std::string indexed(const std::string& property_name, size_t index) const;
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
  /** \brief Associated ParameterHandler
   *
   * If this is well defined, the Properties class can extract parameters
   * directly from the ParameterHandler (only on explicit request via
   * get_property or extract_property).
   */
  ParameterHandler* parameter_handler = nullptr;
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
  /** \brief Extract a property from the parameter handler
   */
  template<typename T>
  bool extract_property(const std::string& property_name, T* property_value);
  /** \brief Account for the different name conventions between ParameterHandler
   * and Property
   */
  std::string parameter_name(const std::string& property_name) const;
};
#endif  // SRC_PROPERTIES_HPP_
