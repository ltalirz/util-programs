#include<iostream>

// Classes inherit all operators of the base class.
// Even if the parameter of the operator is of base class type,
// the operator works for inherited classes (without a defined conversion)


struct A {
    int i;
    A(){ i = 1; }
    A(const A& a2){ i=a2.i; std::cout<< "Copy constructor called\n";}

    // Does not call copy constructor
    const A & operator+=(const A &a2);    
    // Would call copy constructor
    //const A & operator+=(A a2){ i += a2.i; return *this; }

};

struct B : public A {
        int dummy;
        B(){ i=3;}

        // Conversion not needed
        //B(A a)  { i = a.i;}
};

int main(){
    A a = A();
    B b = B();
    b += a;
    a += b;
    std::cout << b.i;
    return 0;
}

const A& A::operator+=(const A &a2){ i += a2.i; return *this; }

