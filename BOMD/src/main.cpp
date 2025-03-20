#include "BOMD.h"

using namespace std;

int main() {
    BOMD bomd("test.gjf", "vel.txt");
    bomd.run(3, 0.5);
}