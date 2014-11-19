/**
 * This file demonstrates the Proxy design pattern.
 * It contains two classes: a "SimpleSet" and a "SimpleElement".
 *
 * When the SimpleSet changes size, all of its elements are copied to new
 * locations in memory and renumbered. But that's OK: since the SimpleElements
 * are proxies, they always access the SimpleSet using the most up-to-date
 * memory address and index.
 *
 * We use dynamic memory allocation here to make the issues obvious; you would
 * almost certainly use an STL container. Our proxy is also slow: it takes O(N)
 * time to access an element, where N is the number of elements. You will fix
 * this in your Graph.
 */

#include <iostream>
#include <string>
#include <cassert>

class SimpleSet {
  // Predeclare the internal struct
  struct internal_element;

 public:
  typedef unsigned size_type;

  /** Constructor
   * Construct a SimpleSet with zero size */
  SimpleSet()
      : elements_(), size_(0), next_uid_(0) {
  }
  /** Destructor */
  ~SimpleSet() {
    delete[] elements_;
  }

  /** The proxy class. This is used to access elements of the SimpleSet class
   *
   * Abstractly, a SimpleElement is a
   *
   */
  class SimpleElement {
   public:
    /** Public Constructor
     * Creates invalid SimpleElement */
    SimpleElement() {
    }
    /** Accessor method to get the text associated with the SimpleElement
     * @pre This SimpleElement is valid and is an element of a SimpleSet
     */
    std::string text(){
      return fetch().text;
    }
    /** Set method which allows for the SimpleElement's text to be changed.
     * This adjusts the value in the SimpleSet as well
     * @pre This SimpleElement is valid and is an element of a SimpleSet
     */
    void set_text(const std::string& text) {
      fetch().text = text;
    }
   private:
    // Pointer back to the SimpleSet container
    SimpleSet* set_;
    // This element's unique identification number
    size_type uid_;
    /** Private Constructor */
    SimpleElement(const SimpleSet* set, size_type uid)
        : set_(const_cast<SimpleSet*>(set)), uid_(uid) {
    }
    /** Helper method to return the appropriate element.
     * This loops over the elements until it finds the element with the
     * correct uid.
     */
    internal_element& fetch() const {
      for (size_type i = 0; i < set_->size(); ++i)
        if (set_->elements_[i].uid == uid_)
          return set_->elements_[i];
      assert(false);
    }
    // Allow SimpleSet to access private variables and methods
    friend class SimpleSet;
  };

  /** Return SimpleSet's size. */
  size_type size() const {
    return size_;
  }
  /** Return a proxy object for element @a i. */
  SimpleElement get_element(size_type i) const {
    assert(i < size());
    return SimpleElement(this, i);
  }
  /** Add a new element at the end.
   * @return A proxy for the new element */
  SimpleElement push_back(const std::string& text) {
    // Create a new elements array
    internal_element* new_elements = new internal_element[size_ + 1];
    // Copy the current elements to a new array
    for (size_type i = 0; i < size_; ++i)
      new_elements[i] = elements_[i];
    // Set the text and uid for the new element
    new_elements[size_].text = text;
    new_elements[size_].uid = next_uid_;
    // Delete the old elements and reassign its value
    delete[] elements_;
    elements_ = new_elements;
    ++size_;
    ++next_uid_;
    // Returns a SimpleElement that points to the new element
    return SimpleElement(this, next_uid_-1);
  }

  /** Remove the element at position @a i, moving later elements down. */
  void remove(size_type i) {
    assert(i < size());
    for (++i; i < size(); ++i)
      elements_[i - 1] = elements_[i];
    --size_;
  }

 private:
  // Internal type for set elements
  struct internal_element {
    std::string text;   // The text held by an element
    size_type uid;      // The unique identifcation for an element
  };

  internal_element* elements_;
  size_type size_;
  size_type next_uid_;

  // Disable copy and assignment of a SimpleSet
  SimpleSet(const SimpleSet&) = delete;
  SimpleSet& operator=(const SimpleSet&) = delete;
};


int main() {
  SimpleSet v;
  SimpleSet::SimpleElement e0 = v.push_back("Hello");
  SimpleSet::SimpleElement e1 = v.push_back("World");
  std::cerr << e0.text() << " " << e1.text() << std::endl;
  // prints "Hello World"

  SimpleSet::SimpleElement e0_copy = v.get_element(0);
  e0.set_text("Goodbye");
  std::cerr << e0.text() << " " << e0_copy.text() << std::endl;
  // prints "Goodbye Goodbye"
  // Since SimpleElement is a proxy, e0 and e0_copy both return
  // the most up-to-date information

  SimpleSet::SimpleElement e2 = v.push_back("Friends");
  v.remove(1);	// This will remove "World" from the simple set
  std::cerr << e0.text() << " " << e2.text() << std::endl;
  // prints "Goodbye Friends": SimpleElement locates its element using a
  // unique number that stays stable even after SimpleSet's internal array
  // is rearranged by the remove

  std::cerr << e1.text() << std::endl;
  // prints Assertion 'false' failed.  This occurs because the SimpleElement e1
  // is no longer a valid element of SimpleSet v: it has been removed.
  // This violates the precondition on text() and causes a failure.
  // Internally, the uid cannot be found in the fetch() method.
}
