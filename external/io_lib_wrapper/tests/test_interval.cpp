// test alignment class
#include "catch.hpp"

#include "interval.hpp"

TEST_CASE( "Test empty interval", "[interval]" ) {
	Interval interval(-1, -1);
	REQUIRE(interval.length() == 0);
}

TEST_CASE( "Test get length", "[interval]" ) {
	Interval interval(10, 11);
	REQUIRE(interval.length() == 2);
}

TEST_CASE( "Test no_overlap", "[interval]" ) {
	Interval interval1(10, 15);
	Interval interval2(16, 20);
	CHECK(interval1.no_overlap(interval2) );
	CHECK(interval2.no_overlap(interval1) );

	Interval interval3(15, 20);
	CHECK_FALSE( interval1.no_overlap(interval3) );
	CHECK_FALSE( interval3.no_overlap(interval1) );
}

TEST_CASE( "Test is_overlapping", "[interval]" ) {
	Interval interval1(10, 15);
	Interval interval2(16, 20);
	CHECK_FALSE(interval1.is_overlapping(interval2) );
	CHECK_FALSE(interval2.is_overlapping(interval1) );

	Interval interval3(15, 20);
	CHECK( interval1.is_overlapping(interval3) );
	CHECK( interval3.is_overlapping(interval1) );
}
