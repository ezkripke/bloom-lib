## bloom-lib
Header-only bloom filter library for C/C++.

Initialize Bloom Filter with total expected number of entries and desired false positive rate:

	BloomFilter bf;
	bloom_init(&bf, 50, 0.05);
	
Add int to filter:

	bloom_update(&bf, 5);
	
Check if int is (probably) in filter
	
	bool b1 = bloom_contains(&bf, 11); // 100% chance that b1 == false
	bool b2 = bloom_contains(&bf, 5);  // 95% chance that b2 == true
