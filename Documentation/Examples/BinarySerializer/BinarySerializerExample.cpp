/**
 * This is an example how to use the BinarySerializer class.
 * The purpose of the BinarySerializer class is to serialize any variable or object to binary raw data.
 * Thanks to this feature one is able to simply dump raw data to a file or send arbitrary data through
 * IPC.
 * This example aims to present how to use the BinarySerializer.
 */
#include <GTMesh/BinarySerializer/BinarySerializer.h>
#include <GTMesh/Debug/Debug.h>
#include <vector>
#include <map>
#include <string>
using namespace std;

namespace ns {

struct Data {
    double d1;
    int i;
};

}
MAKE_ATTRIBUTE_TRAIT(ns::Data, d1, i);

namespace Impl {

void BinarySerialize (BinarySerializer::ByteContainer& cont, const ns::Data& d) {
    BinarySerializer::addNext(cont, 2*d.d1);
    BinarySerializer::addNext(cont, int(0));

}
}

// namespace ns {
// using Impl::BinarySerialize;
// }
template <>
class Interface::Compact<ns::Data, void>: public std::true_type{};

int main()
{
    std::vector<ns::Data> d(10);
    DBGVAR(d);
    BinarySerializer bs;
    bs << d;
    DBGVAR(bs.read<decltype (d)>());
    bs.clear();
    bs << d;


}
