# Implementación del Problema Evolutivo de Navier-Stokes

## Resumen

Se ha implementado la resolución del problema evolutivo (time-dependent) de las ecuaciones de Navier-Stokes en el proyecto MMfem, manteniendo la estructura y estilo del código existente.

## Nuevas Funcionalidades

### 1. Formulaciones en `mmfem/formulations.py`

Se agregaron tres nuevas funciones para manejar problemas dependientes del tiempo:

#### `unsteady_stokes_step()`
Resuelve un paso de tiempo para el problema de Stokes no estacionario usando el método θ:

```python
solution = unsteady_stokes_step(
    mesh, fespace, velocity_prev,
    dirichlet_boundaries, velocity_bc,
    dt=0.01, viscosity=0.01, theta=1.0
)
```

**Parámetros del método θ:**
- `theta = 0.0`: Forward Euler (explícito)
- `theta = 0.5`: Crank-Nicolson (segundo orden)
- `theta = 1.0`: Backward Euler (implícito, incondicionalmente estable) - **default**

#### `unsteady_navier_stokes_semiimplicit_step()`
Realiza un paso de tiempo semi-implícito para Navier-Stokes:
- **Implícito**: derivada temporal + término viscoso (incondicionalmente estable)
- **Explícito**: término convectivo (requiere condición CFL)

```python
solution = unsteady_navier_stokes_semiimplicit_step(
    mesh, fespace, velocity_prev,
    dirichlet_boundaries, velocity_bc,
    dt=0.01, viscosity=0.01,
    convection_form="standard"
)
```

**Ventajas:**
- Rápido (solo un sistema lineal por paso)
- Bueno para números de Reynolds moderados

**Desventajas:**
- Requiere pasos de tiempo pequeños (condición CFL)

#### `unsteady_navier_stokes_implicit_step()`
Realiza un paso de tiempo completamente implícito con linealización:

```python
solution = unsteady_navier_stokes_implicit_step(
    mesh, fespace, velocity_prev, velocity_guess,
    dirichlet_boundaries, velocity_bc,
    dt=0.05, viscosity=0.01,
    method='picard'  # o 'newton'
)
```

**Ventajas:**
- Incondicionalmente estable
- Permite pasos de tiempo más grandes

**Desventajas:**
- Requiere iteraciones no lineales en cada paso de tiempo

### 2. Solvers de Alto Nivel en `mmfem/solvers.py`

#### `time_stepping_semiimplicit()`
Solver completo para evolución temporal con esquema semi-implícito:

```python
solutions, times = time_stepping_semiimplicit(
    mesh, fespace, initial_velocity,
    dirichlet_boundaries, velocity_bc,
    dt=0.01, T_final=1.0,
    viscosity=0.01,
    convection_form="standard",
    save_frequency=10,
    verbose=True,
    callback=my_callback
)
```

#### `time_stepping_implicit()`
Solver completo para evolución temporal con esquema completamente implícito:

```python
solutions, times = time_stepping_implicit(
    mesh, fespace, initial_velocity,
    dirichlet_boundaries, velocity_bc,
    dt=0.05, T_final=1.0,
    viscosity=0.01,
    method='picard',
    tolerance=1e-8,
    max_nonlinear_iter=20,
    save_frequency=5,
    verbose=True
)
```

## Ejemplos de Uso

### Ejemplo 1: Cavity No Estacionario Completo
Archivo: `examples/unsteady_lid_driven_cavity.py`

Simulación completa del problema lid-driven cavity evolutivo con:
- Condición inicial desde solución de Stokes
- Dos opciones de método (semi-implícito o totalmente implícito)
- Cálculo de energía cinética
- Exportación a VTK para visualización en ParaView
- Gráficos de evolución temporal

**Uso:**
```bash
python examples/unsteady_lid_driven_cavity.py
```

### Ejemplo 2: Ejemplo Simple
Archivo: `examples/simple_unsteady_example.py`

Ejemplo minimalista que muestra el uso básico de ambos métodos:
- Flujo con perfil parabólico en la entrada
- Comparación directa de ambos esquemas
- Código simple y educativo

**Uso:**
```bash
python examples/simple_unsteady_example.py
```

## Ecuaciones Implementadas

### Problema de Stokes No Estacionario
```
∂u/∂t - ν Δu + ∇p = 0    en Ω × (0,T)
∇·u = 0                    en Ω × (0,T)
u = u_D                    en ∂Ω × (0,T)
u(x,0) = u_0(x)           en Ω
```

### Problema de Navier-Stokes No Estacionario
```
∂u/∂t - ν Δu + (u·∇)u + ∇p = 0    en Ω × (0,T)
∇·u = 0                             en Ω × (0,T)
u = u_D                             en ∂Ω × (0,T)
u(x,0) = u_0(x)                    en Ω
```

## Esquemas Numéricos

### Semi-Implícito (IMEX)
```
(u^{n+1} - u^n)/Δt - ν Δu^{n+1} + (u^n·∇)u^n + ∇p^{n+1} = 0
```
- Estable para difusión
- Condición CFL: Δt ≤ C·h/||u||

### Totalmente Implícito
```
(u^{n+1} - u^n)/Δt - ν Δu^{n+1} + (u^{n+1}·∇)u^{n+1} + ∇p^{n+1} = 0
```
- Incondicionalmente estable
- Requiere iteraciones Picard o Newton

## Estructura del Código

La implementación sigue la estructura existente del proyecto:

```
MMfem-main/
├── src/mmfem/
│   ├── formulations.py      # ← 3 nuevas funciones agregadas
│   ├── solvers.py            # ← 2 nuevos solvers agregados
│   ├── __init__.py           # ← actualizado con exports
│   ├── mesh.py               # sin cambios
│   └── spaces.py             # sin cambios
├── examples/
│   ├── lid_driven_cavity.py                  # problema estacionario
│   ├── unsteady_lid_driven_cavity.py        # ← NUEVO: problema evolutivo completo
│   └── simple_unsteady_example.py           # ← NUEVO: ejemplo simple
└── tests/
    └── ...
```

## Características Clave

1. **Consistencia con el código existente**: Mismo estilo, convenciones y documentación
2. **Flexibilidad**: Múltiples esquemas temporales disponibles
3. **Callbacks**: Soporte para funciones callback durante la evolución
4. **Verbose output**: Información detallada del progreso
5. **Exportación VTK**: Integración con ParaView para visualización
6. **Documentación completa**: Docstrings detallados con ejemplos

## Uso Básico Paso a Paso

```python
from ngsolve import CF
from mmfem import rectangle_mesh, taylor_hood, time_stepping_semiimplicit
from mmfem.formulations import stokes_problem

# 1. Crear malla
mesh = rectangle_mesh(0, 1, 0, 1, h=0.05)

# 2. Definir espacio de elementos finitos
X = taylor_hood(mesh, "left|right|bottom|top")

# 3. Condiciones de frontera
u_bc = CF((1, 0))

# 4. Condición inicial (solución de Stokes)
stokes_sol = stokes_problem(mesh, X, "top", u_bc, viscosity=0.01)
u0 = stokes_sol.components[0]

# 5. Evolución temporal
solutions, times = time_stepping_semiimplicit(
    mesh=mesh,
    fespace=X,
    initial_velocity=u0,
    dirichlet_boundaries="top",
    velocity_bc=u_bc,
    dt=0.01,
    T_final=1.0,
    viscosity=0.01,
    verbose=True
)

# 6. Usar las soluciones
for i, (sol, t) in enumerate(zip(solutions, times)):
    u, p, lam = sol.components
    # Hacer algo con u, p en el tiempo t
```

## Consideraciones Numéricas

### Elección del Paso de Tiempo

**Semi-implícito:**
- Condición CFL: `dt ≤ C·h/||u||_max`
- Típicamente: `dt ~ 0.01` para `h = 0.05` y velocidad O(1)

**Totalmente implícito:**
- Sin restricción CFL
- Puede usar `dt` 5-10 veces más grande
- Trade-off: iteraciones no lineales

### Elección del Método

| Criterio | Semi-implícito | Totalmente implícito |
|----------|----------------|---------------------|
| Re bajo-moderado | ✓ Preferido | Alternativa |
| Re alto | Usar dt pequeño | ✓ Más robusto |
| Costo computacional | Menor por paso | Mayor por paso |
| Pasos totales | Más pasos | Menos pasos |

## Próximos Pasos Posibles

1. Implementar esquemas de orden superior (BDF2, BDF3)
2. Agregar adaptatividad del paso de tiempo
3. Implementar projection methods (Chorin-Temam)
4. Agregar esquemas de estabilización (SUPG)
5. Paralelización de los solvers temporales

## Referencias

La implementación sigue las formulaciones estándar descritas en:
- Quarteroni & Valli, "Numerical Approximation of Partial Differential Equations" (1994)
- Elman, Silvester & Wathen, "Finite Elements and Fast Iterative Solvers" (2014)
- NGSolve Documentation: https://docu.ngsolve.org/