/* node.cc
 * libxml++ and this file are copyright (C) 2000 by Ari Johnson, and
 * are covered by the GNU Lesser General Public License, which should be
 * included with libxml++ as the file COPYING.
 */

#include <libxml++/nodes/element.h>
#include <libxml++/nodes/node.h>
#include <libxml++/nodes/entityreference.h>
#include <libxml++/nodes/textnode.h>
#include <libxml++/nodes/commentnode.h>
#include <libxml++/nodes/cdatanode.h>
#include <libxml++/nodes/processinginstructionnode.h>
#include <libxml++/exceptions/internal_error.h>
#include <libxml++/document.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/tree.h>

#include <iostream>

namespace xmlpp
{

Node::Node(xmlNode* node)
  : impl_(node)
{
   impl_->_private = this;
}

Node::~Node()
{}

const Element* Node::get_parent() const
{
  return const_cast<Node*>(this)->get_parent();
}

Element* Node::get_parent()
{
  if(!(cobj()->parent && cobj()->parent->type == XML_ELEMENT_NODE))
    return 0;

  Node::create_wrapper(cobj()->parent);
  return static_cast<Element*>(cobj()->parent->_private);
}

const Node* Node::get_next_sibling() const
{
  return const_cast<Node*>(this)->get_next_sibling();
}

Node* Node::get_next_sibling()
{
  if(!cobj()->next)
    return 0;

  Node::create_wrapper(cobj()->next);
  return static_cast<Node*>(cobj()->next->_private);
}

const Node* Node::get_previous_sibling() const
{
  return const_cast<Node*>(this)->get_previous_sibling();
}

Node* Node::get_previous_sibling()
{
  if(!cobj()->prev)
    return 0;

  Node::create_wrapper(cobj()->prev);
  return static_cast<Node*>(cobj()->prev->_private);
}

Node::NodeList Node::get_children(const Glib::ustring& name)
{
   xmlNode* child = impl_->children;
   if(!child)
     return NodeList();

   NodeList children;
   do
   {
      if(name.empty() || name == (const char*)child->name)
      {
        Node::create_wrapper(child);
        children.push_back(reinterpret_cast<Node*>(child->_private));
      }
   }
   while((child = child->next));
   
   return children;
}

const Node::NodeList Node::get_children(const Glib::ustring& name) const
{
  return const_cast<Node*>(this)->get_children(name);
}

Element* Node::add_child(const Glib::ustring& name,
                         const Glib::ustring& ns_prefix)
{
  _xmlNode* child = create_new_child_node(name, ns_prefix);
  if(!child)
    return 0;

  _xmlNode* node = xmlAddChild(impl_, child);
  if(!node)
    return 0;
 
  Node::create_wrapper(node);
  return static_cast<Element*>(node->_private);
}

Element* Node::add_child(xmlpp::Node* previous_sibling, 
                         const Glib::ustring& name,
                         const Glib::ustring& ns_prefix)
{
  if(!previous_sibling)
    return 0;

  _xmlNode* child = create_new_child_node(name, ns_prefix);
  if(!child)
    return 0;

  _xmlNode* node = xmlAddNextSibling(previous_sibling->cobj(), child);
  if(!node)
    return 0;

  Node::create_wrapper(node);
  return static_cast<Element*>(node->_private);
}

Element* Node::add_child_before(xmlpp::Node* next_sibling, 
                         const Glib::ustring& name,
                         const Glib::ustring& ns_prefix)
{
  if(!next_sibling)
    return 0;

  _xmlNode* child = create_new_child_node(name, ns_prefix);
  if(!child)
    return 0;

  _xmlNode* node = xmlAddPrevSibling(next_sibling->cobj(), child);
  if(!node)
    return 0;

  Node::create_wrapper(node);
  return static_cast<Element*>(node->_private);
}

_xmlNode* Node::create_new_child_node(const Glib::ustring& name, const Glib::ustring& ns_prefix)
{
   xmlNs* ns = 0;

   if(impl_->type != XML_ELEMENT_NODE)
   {
      #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
      throw internal_error("You can only add child nodes to element nodes");
      #else
      return 0;
      #endif //LIBXMLCPP_EXCEPTIONS_ENABLED
   }

   //Ignore the namespace if none was specified:
   if(!ns_prefix.empty())
   {
     //Use the existing namespace if one exists:
     ns = xmlSearchNs(impl_->doc, impl_, (const xmlChar*)ns_prefix.c_str());
     if (!ns)
     {
       #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
       throw exception("The namespace prefix (" + ns_prefix + ") has not been declared.");
       #else
       return 0;
       #endif //LIBXMLCPP_EXCEPTIONS_ENABLED
     }
   }

   return xmlNewNode(ns, (const xmlChar*)name.c_str());
}


void Node::remove_child(Node* node)
{
  //TODO: Allow a node to be removed without deleting it, to allow it to be moved?
  //This would require a more complex memory management API.
  
  xmlNode* cnode = node->cobj();
  Node::free_wrappers(cnode); //This delete the C++ node (not this) itself.
  xmlUnlinkNode(cnode);
  xmlFreeNode(cnode);
}

Node* Node::import_node(const Node* node, bool recursive)
{
  //Create the node, by copying:
  xmlNode* imported_node = xmlDocCopyNode(const_cast<xmlNode*>(node->cobj()), impl_->doc, recursive);
  if (!imported_node)
  {
    #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
    throw exception("Unable to import node");
    #else
    return 0;
    #endif //LIBXMLCPP_EXCEPTIONS_ENABLED
  }

  //Add the node:
  xmlNode* added_node = xmlAddChild(this->cobj(),imported_node);
  if (!added_node)
  {
    Node::free_wrappers(imported_node);
    xmlFreeNode(imported_node);

    #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
    throw exception("Unable to add imported node to current node");
    #else
    return 0;
    #endif //LIBXMLCPP_EXCEPTIONS_ENABLED
  }

  Node::create_wrapper(imported_node);
  return static_cast<Node*>(imported_node->_private);
}

Glib::ustring Node::get_name() const
{
  return impl_->name ? (const char*)impl_->name : "";
}

void Node::set_name(const Glib::ustring& name)
{
  xmlNodeSetName( impl_, (const xmlChar *)name.c_str() );
}

int Node::get_line() const
{
   return XML_GET_LINE(impl_);
}


xmlNode* Node::cobj()
{
  return impl_;
}

const xmlNode* Node::cobj() const
{
  return impl_;
}

Glib::ustring Node::get_path() const
{
  xmlChar* path = xmlGetNodePath(impl_);
  Glib::ustring retn = path ? (char*)path : "";
  xmlFree(path);
  return retn;
}

static NodeSet find_impl(xmlXPathContext* ctxt, const Glib::ustring& xpath)
{
  xmlXPathObject* result = xmlXPathEval((const xmlChar*)xpath.c_str(), ctxt);

  if(!result)
  {
    xmlXPathFreeContext(ctxt);

    #ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
    throw exception("Invalid XPath: " + xpath);
    #else
    return NodeSet();
    #endif //LIBXMLCPP_EXCEPTIONS_ENABLED
  }

  if(result->type != XPATH_NODESET)
  {
    xmlXPathFreeObject(result);
    xmlXPathFreeContext(ctxt);

    #ifdef LIBXMLCPP_EXCEPTIONS_ENABLE
    throw internal_error("Only nodeset result types are supported.");
    #else
    return NodeSet();
    #endif //LIBXMLCPP_EXCEPTIONS_ENABLED
  }

  xmlNodeSet* nodeset = result->nodesetval;
  NodeSet nodes;
  if( nodeset && !xmlXPathNodeSetIsEmpty(nodeset))
  {
    const int count = xmlXPathNodeSetGetLength(nodeset);
    nodes.reserve(count);
    for (int i = 0; i != count; ++i)
    {
      xmlNode* cnode = xmlXPathNodeSetItem(nodeset, i);
      if(cnode->type == XML_NAMESPACE_DECL)
      {
        //In this case we would cast it to a xmlNs*,
        //but this C++ method only returns Nodes.
        std::cerr << "Node::find_impl: ignoring an xmlNs object." << std::endl;
        continue;
      }
      
      //TODO: Check for other cnode->type values?
  
      Node::create_wrapper(cnode);
      Node* cppNode = static_cast<Node*>(cnode->_private);
      nodes.push_back(cppNode);
    }
  }
  else
  {
    // return empty set
  }

  xmlXPathFreeObject(result);
  xmlXPathFreeContext(ctxt);

  return nodes;
}

NodeSet Node::find(const Glib::ustring& xpath) const
{
  xmlXPathContext* ctxt = xmlXPathNewContext(impl_->doc);
  ctxt->node = impl_;
  
  return find_impl(ctxt, xpath);
}

NodeSet Node::find(const Glib::ustring& xpath,
		   const PrefixNsMap& namespaces) const
{
  xmlXPathContext* ctxt = xmlXPathNewContext(impl_->doc);
  ctxt->node = impl_;

  for (PrefixNsMap::const_iterator it=namespaces.begin();
       it != namespaces.end(); it++)
    xmlXPathRegisterNs(ctxt,
		       reinterpret_cast<const xmlChar*>(it->first.c_str()),
		       reinterpret_cast<const xmlChar*>(it->second.c_str()));

  return find_impl(ctxt, xpath);
}

Glib::ustring Node::get_namespace_prefix() const
{
  if(impl_->type == XML_DOCUMENT_NODE)
  {
    //impl_ is actually of type xmlDoc, instead of just xmlNode.
    //libxml does not always use GObject-style inheritance, so xmlDoc does not have all the same struct fields as xmlNode.
    //Therefore, a call to impl_->ns would be invalid.
    //This can be an issue when calling this method on a Node returned by Node::find().
    //See the TODO comment on Document, suggesting that Document should derived from Node.

    return Glib::ustring();
  }

  if(impl_ && impl_->ns && impl_->ns->prefix)
    return (char*)impl_->ns->prefix;
  else
    return Glib::ustring();
}

Glib::ustring Node::get_namespace_uri() const
{
  if(impl_->type == XML_DOCUMENT_NODE)
  {
    //impl_ is actually of type xmlDoc, instead of just xmlNode.
    //libxml does not always use GObject-style inheritance, so xmlDoc does not have all the same struct fields as xmlNode.
    //Therefore, a call to impl_->ns would be invalid.
    //This can be an issue when calling this method on a Node returned by Node::find().
    //See the TODO comment on Document, suggesting that Document should derived from Node.

    return Glib::ustring();
  }

  if(impl_ && impl_->ns && impl_->ns->href)
    return (char*)impl_->ns->href;
  else
    return Glib::ustring();
}

void Node::set_namespace(const Glib::ustring& ns_prefix)
{
  //Look for the existing namespace to use:
  xmlNs* ns = xmlSearchNs( cobj()->doc, cobj(), (xmlChar*)(ns_prefix.empty() ? 0 : ns_prefix.c_str()) );
  if(ns)
  {
      //Use it for this element:
      xmlSetNs(cobj(), ns);
  }
  else
  {
    #ifdef LIBXMLCPP_EXCEPTIONS_ENABLE
    throw exception("The namespace (" + ns_prefix + ") has not been declared.");
    #endif //LIBXMLCPP_EXCEPTIONS_ENABLE
  }
}

void Node::create_wrapper(xmlNode* node)
{
  if(node->_private)
  {
	  //Node already wrapped, skip
	  return;
  }

  switch (node->type)
  {
    case XML_ELEMENT_NODE:
    {
      node->_private = new xmlpp::Element(node);
      break;
    }
    case XML_ATTRIBUTE_NODE:
    {
      node->_private = new xmlpp::Attribute(node);
      break;
    }
    case XML_TEXT_NODE:
    {
      node->_private = new xmlpp::TextNode(node);
      break;
    }
    case XML_COMMENT_NODE:
    {
      node->_private = new xmlpp::CommentNode(node);
      break;
    }
    case XML_CDATA_SECTION_NODE:
    {
      node->_private = new xmlpp::CdataNode(node);
      break;
    }
    case XML_PI_NODE:
    {
      node->_private = new xmlpp::ProcessingInstructionNode(node);
      break;
    }
    case XML_DTD_NODE:
    {
      node->_private = new xmlpp::Dtd(reinterpret_cast<xmlDtd*>(node));
      break;
    }
    //case XML_ENTITY_NODE:
    //{
    //  assert(0 && "Warning: XML_ENTITY_NODE not implemented");
    //  //node->_private = new xmlpp::ProcessingInstructionNode(node);
    //  break;
    //}
    case XML_ENTITY_REF_NODE:
    {
      node->_private = new xmlpp::EntityReference(node);
      break;
    }
    case XML_DOCUMENT_NODE:
    {
      // do nothing. For Documents it's the wrapper that is the owner.
      break;
    }
    default:
    {
      // good default for release versions
      node->_private = new xmlpp::Node(node);
      std::cerr << G_STRFUNC << "Warning: new node of unknown type created: " << node->type << std::endl;
      break;
    }
  }
}

void Node::free_wrappers(xmlNode* node)
{
  if(!node)
    return;
    
  //Walk the children list
  for(xmlNode* child=node->children; child; child=child->next)
    free_wrappers(child);

  //Delete the local one
  switch(node->type)
  {
    //Node types that have no properties
    case XML_DTD_NODE:
      delete static_cast<Dtd*>(node->_private);
      node->_private = 0;
      return;
    case XML_ATTRIBUTE_NODE:
    case XML_ELEMENT_DECL:
    case XML_ATTRIBUTE_DECL:
    case XML_ENTITY_DECL:
      delete static_cast<Node*>(node->_private);
      node->_private = 0;
      return;
    case XML_DOCUMENT_NODE:
      //Do not free now. The Document is usually the one who owns the caller.
      return;
    default:
      delete static_cast<Node*>(node->_private);
      node->_private = 0;
      break;
  }

  //Walk the attributes list.
  //Note that some "derived" struct have a different layout, so 
  //_xmlNode::properties would be a nonsense value, leading to crashes,
  //(and shown as valgrind warnings), so we return above, to avoid 
  //checking it here.
  for(xmlAttr* attr = node->properties; attr; attr = attr->next)
    free_wrappers(reinterpret_cast<xmlNode*>(attr));
}


} //namespace xmlpp
