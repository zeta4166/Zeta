# Hi. I'm Zeta.

Synchrotron science and materials engineering.

![BLAK'AT](assets/blakat.gif)

More about me...

```python
from dataclasses import dataclass, field

@dataclass
class Zeta:
    job: str = "Engineer"
    favorite_fields: list = field(default_factory=lambda: [
        "X-ray Diffraction", "Spectroscopy", "Condensed Matter", "Materials Science"
    ])
    main_tools: list = field(default_factory=lambda: [
        "Python", "NumPy", "MATLAB", "LaTeX"
    ])
    favorite_language: list = field(default_factory=lambda: [
        "Python", "C", "JavaScript"
    ])

    def state(self, t):
        # Sit and drink pennyroyal tea
        return self

# I was an ordinary person who studied hard. There are no miracle people.
```

## What I'm working on

- Numerical simulations of physical systems
- Engineering projects with a healthy dose of physics
- Learning something new about the universe every day

## Fun fact

> "If you can't explain it simply, you don't understand it well enough."

I take that as a challenge.
