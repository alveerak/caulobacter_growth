---
layout: page
title: Analysis
permalink: interactive_a
sidebar: true
interactive: interactive_1.html
---
---



{% for entry in site.data.analysis %}

{% if entry[0] != 'title' %}
{% if entry[0] != 'authors' %}
## {{entry[0]}}
{{entry[1]}}

{% endif %}
{% endif %}

{% if entry[0] == 'Choosing and Performing Segmentation on a Caulobacter crescentus cell' %}
<img src="{{ site.baseurl }}/assets/img/fig1.gif" width="500" height="400">
<img src="{{ site.baseurl }}/assets/img/fig2.gif" width="500" height="400">

{% endif %}

{% if entry[0] == 'Graphing Area over time for both dividing bacteria' %}

<img src="{{ site.baseurl }}/assets/img/graph1.png" width="500" height="400">
<img src="{{ site.baseurl }}/assets/img/graph2.png" width="900" height="350">
{% endif %}

{% if entry[0] == 'Color Code Growth Events' %}
<img src="{{ site.baseurl }}/assets/img/graph3.png" width="500" height="400">
<img src="{{ site.baseurl }}/assets/img/graph4.png" width="1000" height="300">

{% endif %}

{% if entry[0] == 'ECDFs and Analysis of Time Difference' %}

<img src="{{ site.baseurl }}/assets/img/graph5.png" width="500" height="400">
<img src="{{ site.baseurl }}/assets/img/graph6.png" width="500" height="400">

{% endif %}

{% if entry[0] == 'Acknowledgments' %}

{% endif %}
{% endfor %}  
