## MATH 3610 (Introduction to Mathematical Modelling)

This GitHub repo contains material related to the University of Manitoba course MATH 3610, Introduction to Mathematical Modelling and reflects the organisation of the course when I ([Julien Arino](https://julien-arino.github.io/)) teach it.

Here, you will find the publicly available information: code, slides and some general information. All other information (syllabus, assignments, marks) is available through UM Learn.

### This GitHub repo

On the [GitHub version](https://github.com/julien-arino/math-3610-introduction-to-math-modelling/) of the page, you have access to all the files. You can also download the entire repository by clicking the buttons on the left. (You can also of course clone this repo, but you will need to do that from the GitHub version of the site.)

Feel free to use the material in these slides or in the folders. If you find this useful, I will be happy to know.

### Slides

Note that, at present, the slides are work in progress and are only checked up to the point covered in lectures.

<ul>
{% for file in site.static_files %}
  {% if file.path contains 'SLIDES' %}
    {% if file.path contains 'math-3610' %}
      {% if file.path contains 'pdf' %}
        <li><a href="https://julien-arino.github.io/math-3610-introduction-to-math-modelling/SLIDES/{{ file.basename }}.pdf">{{ file.basename }}</a></li>
      {% endif %}
    {% endif %}
  {% endif %}
{% endfor %}
</ul>

### Additional slides and videos

The following slides and videos are part of a set of "vignettes" about R (see [here](https://julien-arino.github.io/R-for-modellers/)). They are required reading/watching for the course as they explain some mechanisms that your R assignments will need to use.

- [Installing and loading packages](https://julien-arino.github.io/R-for-modellers/SLIDES/vignette-03-installing-using-packages.html) ([video](https://youtu.be/WPYHU2G7U7Q))