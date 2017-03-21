SUBDIRS = gctp src 

RECURSIVE_TARGETS = all-recursive clean-recursive install-recursive copy-makefile-recursive

all: all-recursive

install: install-recursive

clean: clean-recursive

copy-makefile: copy-makefile-recursive

$(RECURSIVE_TARGETS):
	@target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir ..."; \
	  (cd $$subdir && $(MAKE) $$target); \
	done;

