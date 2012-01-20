#include <iostream>
#include <vector>

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
    pushBack();
    return 0;
}
