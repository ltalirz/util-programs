#include <iostream>
#include <vector>
#include <iterator>
#include <sstream>



// Out of range iterator alone does not give segmentation fault.
void testIterator(){
    std::vector<int> v(2,10);
    std::vector<int>::iterator it = v.begin();
    ++it;++it;++it;++it;
}




// I experienced problems passing stringstreams by reference
// (after some correct data, there was bogus).
// I am not able to reproduce them here. 
void resizeReference(std::stringstream &s){
    // Fill with 1E7 * 4 Bytes = 40 MByte of data
    int n = 10E6;
    for(int i = 0; i<n; ++i){
        s << "a";
    }
}
void passReference(){
    std::stringstream s;
    resizeReference(s);
    std::string str = s.str();
    std::cout << str.size();
}


// Operator << is not overloaded to cout vector
void coutVector(){
    std::vector<int> v; v.push_back(1);
    //std::cout << v;
}


class DefaultConstr {
    public:
    std::vector<int> v;
};

// Testing the default initialization of a vector
// Calling size() on an empty vector is not a problem
// apparently.
void testDefaultConstr(){
    DefaultConstr d = DefaultConstr();
    std::cout << d.v.size();
}



// Returning a const reference to a class member is problematic, since the
// object may be destroyed before the reference.
// Without new/delete, however, I think one cannot break it.
// The alternative is to define friend classes/functions 
class Test {
    public:
        double d;
        Test(): d(5) {};
        void setD(double n){ d=n;}
        const double & getD() {return d;}
};

void constReference() {
    Test *t = new Test();
    const double &d1 =  t->getD();
    double d2 = t->getD();
    t->setD(4);
    std::cout << d1;
    std::cout << d2;
    delete t;
    std::cout << d1;  // Writes 0

}

// You may pass iterators by reference
void iteratorReference(std::string::iterator &it) {
    ++it;
}

// Double pushback is not allowed.
void pushBack() {
    std::vector<int> numbers;
//  numbers.push_back(3).push_back(4);

}

// Assignment calls the copy constructor
void copyVector() {
    std::vector<int> numbers;
    numbers.push_back(1231);
    std::vector<int> newNumbers = numbers;
    numbers[0]=120;
    std::cout << newNumbers[0];
}

// Reserve leaves elements untouched
void reserveVector() {
    std::vector<int> numbers;
    numbers.push_back(1231);
    numbers.reserve(3);
    std::cout << numbers[0];
    numbers[1]=3;
}

// operator== on vector
void compareVectors() {
    std::vector<int> v1, v2;
    v1.push_back(1);
    v2.push_back(1);
    std::cout << (v1 == v2);
}

int main() {
    //reserveVector();
    //pushBack();
    //std::string s = "HI";
    //std::string::iterator it = s.begin();
    //iteratorReference(it);
    //std::cout << *it;
    //constReference();
    //compareVectors();
    //testDefaultConstr();
    //coutVector();
    //passReference();
    //doVector();
    testIterator();
    return 0;
}
