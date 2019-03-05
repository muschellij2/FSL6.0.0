#include <iostream>
#include <stdexcept>
#include <glibmm/ustring.h>
#include <cstdlib>
#include <libxml++/libxml++.h>

using namespace xmlpp;
using namespace std;

int main (int /* argc */, char** /* argv */)
{
  // Set the global C and C++ locale to the user-configured locale,
  // so we can use std::cout with UTF-8, via Glib::ustring, without exceptions.
  std::locale::global(std::locale(""));

  #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
  try
  {
  #endif //LIBXMLCPP_EXCEPTIONS_ENABLED 
    DomParser example1("example1.xml");
    DomParser example2("example2.xml");
    
    Document *doc1 = example1.get_document();
    Document *doc2 = example2.get_document();
    
    Element *root1 = doc1->get_root_node();
    Element *root2 = doc2->get_root_node();

    // find the first "child" element in example2
    Node::NodeList child_list = root2->get_children("child");
    Node *node_to_add = child_list.front();

    // import the node under the root element (recursive is default)
    root1->import_node(node_to_add);
    
    // print out the new doc1
    string doc1_string = doc1->write_to_string_formatted();
    cout << doc1_string;
    return EXIT_SUCCESS;
  #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
  }
  catch (std::exception& ex)
  {
    cerr << "Caught exception " << ex.what() << endl;
    return EXIT_FAILURE;
  }
  #endif //LIBXMLCPP_EXCEPTIONS_ENABLED 
}
