# Research Website for Caulobacter Growth&Division Analysis
Cloning this repository will allow you to locally load up a website with resources and analysis regarding 
Caulobacter crescentus cell division.

It was made by modifying a bare-bones template website that allowed for publishing of our work as
an interactive website. It's tailored for easy inclusion of data files, plots,
and interactive figures. 

This repository is designed to be hosted as a website on GitHub using the
[GitHub Pages]() hosting service. This service is based on the Ruby framework
[Jekyll]() which is tailored for blog management. I've made several
modifications to the structure to be amenable for scientific publications. 

## Software requirements 
To deploy locally, you must have a [Ruby development environment]() installed as
well as [Jekyll](). You can install Jekyll via 

```
gem install jekyll bundler
```

Once you have Jekyll installed, you can install of the Ruby requirements for
this website by running the following in the command line from the template directory:

```
bundle install
```

The build and preview the website locally, execute the following:
```
bundle exec jekyll serve --watch
```

This will build the website and serve it up at the address:
[http://127.0.0.1:4000](http://127.0.0.1:4000).


## License
This template is a heavily modified version of [Flexible
Jekyll](https://artemsheludko.github.io/flexible-jekyll/) under a GNU General
Public License Version 3.0. This template is provided with the same license.
All writing, logo, and other creative works provided with this template are
issued with a [CC-BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).


```
Copyright (C) 2019  Griffin Chure 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```
