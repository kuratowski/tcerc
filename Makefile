MAKE = make

DIRS = $(addprefix build/,src lib tests examples)

all check demo: $(DIRS) build/Makefile
	cd build && $(MAKE) $@
$(DIRS):
	mkdir -p $@
build/Makefile: Makefile.build
	cp $< $@

clean:
	rm -rf build
