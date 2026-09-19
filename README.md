# Hi. I'm Zeta.

Engineer who thinks in equations.

![BLAK'AT](https://makeagif.com/gif/crystal-castles-kept-U5flek)

More about me...

```python
from dataclasses import dataclass, field

@dataclass
class Zeta:
    job: str = "Engineer"
    pronouns: str = "He/Him"
    passion: str = "Physics"
    favorite_fields: list = field(default_factory=lambda: [
        "Classical Mechanics", "Electromagnetism", "Quantum Mechanics"
    ])
    main_tools: list = field(default_factory=lambda: [
        "Python", "NumPy", "MATLAB", "LaTeX"
    ])
    favorite_language: list = field(default_factory=lambda: [
        "Python", "C", "JavaScript"
    ])
    favorite_equation: str = "∇·E = ρ/ε₀"

    def state(self, t):
        # Every system evolves; I just try to write down the Hamiltonian.
        return self

# As a kid, I was taught that the universe runs on math. I never got over it.
```

## What I'm working on

- Numerical simulations of physical systems
- Engineering projects with a healthy dose of physics
- Learning something new about the universe every day

## Fun fact

> "If you can't explain it simply, you don't understand it well enough."

I take that as a challenge.
