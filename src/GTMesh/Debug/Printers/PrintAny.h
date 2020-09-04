#ifndef PRINTANY_H
#define PRINTANY_H

#include <iomanip>
struct PrintAny {

    template <typename T>
    static void print(std::ostream& ost, const T& var, ...) {
        ost << typeid (T).name() << "(";
        ost << std::hex;
        for (char* byte = (char*)&var; byte < (char*)&var + sizeof (T); byte++) {
            ost << std::setfill('0') << std::setw(2) << ((int)(*byte) & 0xff);
            if (byte < (char*)&var + sizeof (T) - 1) {
                ost << ' ';
            }
        }
        ost << std::dec;
        ost << ")";
    }

    template <typename T>
    static void print(const T& var, ...) {
        for (char* byte = (char*)&var; byte < (char*)&var + sizeof (T); byte++) {
            printf("%02X",((int)(*byte) & 0xff));
            if (byte < (char*)&var + sizeof (T) - 1) {
                printf(" ");
            }
        }
    }

};

#endif // PRINTANY_H
