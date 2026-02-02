# Guía de Visualización en ParaView

## 🎬 Cómo ver animaciones de tus simulaciones

Después de ejecutar una simulación evolutiva con MMfem, obtendrás archivos VTK que puedes visualizar en ParaView.

## 📦 Archivos Generados

Cuando ejecutas una simulación, se generan:

1. **Archivo PVD** (e.g., `cavity_animation.pvd`)
   - Archivo de colección que contiene toda la serie temporal
   - **Este es el archivo que debes abrir en ParaView**

2. **Carpeta con archivos VTU** (e.g., `vtk_output/`)
   - Contiene un archivo `.vtu` para cada paso de tiempo
   - ParaView los lee automáticamente a través del archivo PVD

## 🚀 Pasos para Visualizar en ParaView

### 1. Abrir el archivo PVD

```bash
# Desde la línea de comandos
paraview cavity_animation.pvd

# O desde ParaView:
# File → Open → Seleccionar cavity_animation.pvd
```

### 2. Cargar los datos

1. En el panel "Properties", haz click en **"Apply"**
2. Los datos se cargarán y verás la primera solución

### 3. Reproducir la animación

**Controles de animación** (parte superior de ParaView):
- ▶️ **Play**: Reproduce la animación
- ⏸️ **Pause**: Pausa la animación
- ⏮️ **First Frame**: Primer paso de tiempo
- ⏭️ **Last Frame**: Último paso de tiempo
- 🎚️ **Time slider**: Arrastra para navegar manualmente

### 4. Visualizar el campo de velocidad

**Para ver vectores de velocidad:**

1. En el pipeline (izquierda), selecciona el objeto cargado
2. En "Properties", selecciona "velocity" en el menú desplegable
3. Click **"Apply"**
4. Agrega flechas:
   - Filters → Common → **Glyph**
   - Glyph Type: Arrow
   - Scale Mode: vector
   - Scale Factor: 0.1 (ajusta según necesites)
   - Click "Apply"

**Para ver magnitud de velocidad con colores:**

1. En "Coloring", selecciona "velocity" → "Magnitude"
2. Ajusta la escala de colores en "Color Map Editor"

### 5. Visualizar la presión

1. En "Coloring", cambia a "pressure"
2. Ajusta la escala de colores según necesites

### 6. Visualizar líneas de corriente (streamlines)

1. Filters → Common → **Stream Tracer**
2. Seed Type: Line Source o Point Cloud
3. Vectors: velocity
4. Click "Apply"
5. Opcional: Filters → Tube para hacer las líneas más visibles

## 📊 Visualizaciones Recomendadas

### Para Cavity Flow (flujo en cavidad):

```python
# Después de la simulación:
from mmfem import export_time_series

pvd_file = export_time_series(
    mesh=mesh,
    solutions=solutions,
    times=times,
    field_names=["velocity", "pressure"],
    output_name="cavity_flow",
    verbose=True
)
```

**En ParaView:**
1. Abrir `cavity_flow.pvd`
2. Apply
3. Coloring → velocity → Magnitude
4. Add Filter → Glyph (flechas de velocidad)
5. Add Filter → Stream Tracer (líneas de corriente)
6. Play animation

### Para Channel Flow (flujo en canal):

```python
pvd_file = export_time_series(
    mesh=mesh,
    solutions=solutions,
    times=times,
    field_names=["velocity", "pressure"],
    output_name="channel_flow",
    verbose=True
)
```

**En ParaView:**
1. Abrir `channel_flow.pvd`
2. Apply
3. Slice perpendicular al flujo:
   - Filters → Common → **Slice**
   - Normal: [1, 0, 0] (perpendicular al eje x)
4. Coloring → velocity → Magnitude
5. Ver evolución de perfiles de velocidad

## 🎨 Tips de Visualización

### Mejorar la apariencia

1. **Superficie suave:**
   - Filters → Alphabetical → **Extract Surface**
   - Representation: Surface With Edges

2. **Escala de colores:**
   - Click en "Edit Color Map" (icono de paleta)
   - Choose Preset: "Cool to Warm", "Viridis", etc.
   - Enable "Use Log Scale" si hay grandes variaciones

3. **Anotaciones:**
   - Sources → **Text**
   - Escribe título o tiempo actual
   - Usa "Time Annotation" para mostrar tiempo dinámicamente

### Exportar animaciones

**Como imágenes:**
1. File → Save Animation
2. Selecciona formato (PNG, JPEG)
3. Configura framerate y resolución
4. Save

**Como video:**
1. File → Save Animation
2. Selecciona formato (AVI, OGV)
3. Configura framerate
4. Save

## 🐛 Problemas Comunes

### "No se ve la animación, solo un frame"

❌ **Problema:** Abriste un archivo `.vtu` individual
✅ **Solución:** Abre el archivo `.pvd` en su lugar

### "Los vectores son muy pequeños/grandes"

✅ **Solución:** En Glyph, ajusta "Scale Factor"

### "No veo las flechas de velocidad"

✅ **Solución:** 
1. Asegúrate de aplicar el filtro Glyph
2. En Glyph, verifica que "Vectors" esté en "velocity"
3. Aumenta "Scale Factor"

### "La animación va muy rápido/lento"

✅ **Solución:**
- View → Animation View
- Ajusta "Duration" (segundos totales)
- O ajusta framerate en preferencias

## 💡 Ejemplo Completo de Workflow

```python
# 1. Ejecutar simulación
from mmfem import (
    rectangle_mesh, taylor_hood, stokes_problem,
    time_stepping_semiimplicit, export_time_series
)
from ngsolve import CF

mesh = rectangle_mesh(0, 1, 0, 1, h=0.05)
X = taylor_hood(mesh, "left|right|bottom|top")
u_bc = CF((1, 0))

stokes_sol = stokes_problem(mesh, X, "top", u_bc, viscosity=0.01)
u0 = stokes_sol.components[0]

solutions, times = time_stepping_semiimplicit(
    mesh, X, u0, "top", u_bc,
    dt=0.01, T_final=1.0, viscosity=0.01,
    save_frequency=5
)

# 2. Exportar para ParaView
pvd_file = export_time_series(
    mesh, solutions, times,
    field_names=["velocity", "pressure"],
    output_name="mi_simulacion"
)

print(f"Abre con: paraview {pvd_file}")
```

**En ParaView:**
```
1. paraview mi_simulacion.pvd
2. Apply
3. Coloring → velocity → Magnitude
4. Filters → Glyph (vectores)
5. Play ▶️
```

## 📚 Recursos Adicionales

- [ParaView Tutorial](https://www.paraview.org/Wiki/The_ParaView_Tutorial)
- [ParaView Guide](https://docs.paraview.org/en/latest/)
- [NGSolve VTK Output](https://docu.ngsolve.org/latest/i-tutorials/unit-6.1.2-vtkout/vtkout.html)
