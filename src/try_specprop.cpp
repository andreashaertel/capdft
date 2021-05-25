
#include <iostream>
#include <vector>
#include <typeinfo>



class data {
  public: 
    const std::type_info *my_type;
    virtual ~data() {};
  //get_data
};

template <class T> 
class my_data : public data {
  public: 
    T value;
    my_data(T new_value) {
      value = new_value;
      this->my_type = &(typeid(T));
    };
};




class properties {
  std::vector<data*> myvector;

  public: 
    properties(void) {
      myvector.clear();
    };
    template <class T> void add_property(T value) {
      my_data<T> *new_value = new my_data<T>(value);
      myvector.push_back(new_value);
    };
    template <class T> T get_property(void) {
      // check data type
      if (typeid(T) != *(myvector.back()->my_type)) { std::cout << "moep" << std::endl; return 0; }
      return (dynamic_cast<my_data<T>*>(myvector.back()))->value;
    };
};



int main(int argc, char** argv) {
  my_data<int> *value_a = new my_data<int>(5);
  my_data<int> *value_b = new my_data<int>(4);
  my_data<double> *value_c = new my_data<double>(4.5);

  //myvector.push_back(value_a);
  //myvector.push_back(value_b);
  //myvector.push_back(value_c);

  data *val = value_c;
  
  //my_data<double> *value_d = dynamic_cast<my_data<double>*>(myvector.back());
  my_data<double> *value_d = dynamic_cast<my_data<double>*>(val);

  std::cout << " C: " << value_d->value << std::endl;

  properties *props = new properties();

  props->add_property<double>(3.1);
  std::cout << " Get: " << props->get_property<int>() << std::endl;
  std::cout << " Get: " << props->get_property<double>() << std::endl;

  delete value_c;
  delete value_b;
  delete value_a;
  return 0;
}

