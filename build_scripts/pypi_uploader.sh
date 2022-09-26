cd ..
python setup.py sdist bdist_wheel
twine upload dist/* -r splotch
rm -rf build/
rm -rf dist/
rm -rf splotch.egg-info/
