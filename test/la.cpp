#include "la.hpp"
#include "types.hpp"
#include <vector>
#include <iterator>

// Returning a const reference to a class member is problematic, since the
// object may be destroyed before the reference.
// Without new/delete, however, I think one cannot break it.
// The alternative is to define friend classes/functions 
void testCell(){
    using namespace la;
    Direction d = Direction();
    std::vector<types::Real> v;
    v.push_back(1.0);
    v.push_back(1.0);
    v.push_back(1.0);

    d.incrementVector =v;
    d.incrementCount =5;
    std::vector<Direction> directions;
    directions.push_back(d);
    Cell c = Cell(directions);
}
int main() {
    testCell();
    return 0;
}
