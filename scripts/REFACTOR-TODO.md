TODO
====

1. Teach dynein.run to take paper_params (or other params hypothetically) as an option.

2. Teach dynein.run to use rq?

3. Rewrite generate-paper-*.py to use dynein.run, not os.system('rq run ...')

4. Write single script to regenerate *all* paper data.

5. Test single script by deleting all data and then regenerating it.

6. Do not rename files for paper.

7. Teach plotting scripts to read paper_params module and figure out
   the proper filenames themselves.
