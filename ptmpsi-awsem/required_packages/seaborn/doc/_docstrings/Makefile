rst_files := $(patsubst %.ipynb,../docstrings/%.rst,$(wildcard *.ipynb))
export MPLBACKEND := module://matplotlib_inline.backend_inline

docstrings: ${rst_files}

../docstrings/%.rst: %.ipynb
	../tools/nb_to_doc.py $*.ipynb ../docstrings
	@cp -r ../docstrings/$*_files ../generated/
	@if [ -f ../generated/seaborn.$*.rst ]; then \
	    touch ../generated/seaborn.$*.rst; \
	 fi

clean:
	rm -rf ../docstrings
