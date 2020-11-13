all: 
	sh/build.sh

test: all
	build/libsemigroups/test_bmat8_enum "[quick]"

full-3: 
	@sh/full-run.sh 3

full-4: 
	@sh/full-run.sh 4

full-5: 
	@sh/full-run.sh 5

full-6: 
	@sh/full-run.sh 6

full-7: 
	@sh/full-run.sh 7 8

full-8: 
	@sh/full-run.sh 8 

ref-3: 
	@sh/ref-run.sh 3

ref-4: 
	@sh/ref-run.sh 4

ref-5: 
	@sh/ref-run.sh 5

ref-6: 
	@sh/ref-run.sh 6

clean:  
	rm -rf build
