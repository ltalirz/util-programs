#include <iostream>
#include <vector>
#include <iterator>

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

int main() {
    //reserveVector();
    //pushBack();
    //std::string s = "HI";
    //std::string::iterator it = s.begin();
    //iteratorReference(it);
    //std::cout << *it;
    constReference();
    return 0;
}
