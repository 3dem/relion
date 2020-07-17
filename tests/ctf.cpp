#include <catch2/catch.hpp>
#include "src/ctf.h"

//Actually test the getCTF function. You may wish to test the CTF constructor and setters/getters separately.
TEST_CASE( "Test getCTF", "[ctf]" ) {
  CTF ctf;
  ctf.setValues(10000.0, 12000.0, 90.0, 300.0, 2.7, 0.1, 0.0, 1.0, 0.0);
  float val = ctf.getCTF(10.0, 10.0);
  REQUIRE(val == Approx(0.59154));
}
