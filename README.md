	BloomFilter bf;
	bloom_init(&bf, 50, 0.05);
	for (int i = 0; i < 100; i+=2) {
		bloom_update(&bf, i);
	}
	bool b1 = bloom_contains(&bf, 17); // 100% chance that b1 == false
	bool b2 = bloom_contains(&bf, 8); // 95% chance that b2 == true
