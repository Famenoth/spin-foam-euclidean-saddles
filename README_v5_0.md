# Spin-Foam Euclidean Saddles - GUI v5.0 ✨

Interface gráfica completa para análise numérica do artigo:
**"Causally Confined Euclidean Saddles in Spin-Foam Quantum Gravity"**

## 🎉 NEW IN v5.0: 3D Light Cone Visualization!

Agora com visualização interativa completa do cone de luz de Minkowski, vetores normais, e restrições de fechamento!

## Instalação

### Requisitos

```bash
pip install numpy scipy matplotlib PyQt5
```

### Versões testadas
- Python 3.8+
- NumPy 1.21+
- SciPy 1.7+
- Matplotlib 3.4+
- PyQt5 5.15+

## Executando a Aplicação

```bash
python spin_foam_gui_v5_0.py
```

## 🆕 Funcionalidades v5.0

### 🔮 **Tab: 3D Visualization**

Visualização interativa 3D completa do espaço de Minkowski:

#### **Elementos Visuais:**

**1. Cone de Luz:**
- Cone futuro (cyan/azul claro)
- Cone passado (magenta)
- Bordas geométricas

**2. Hiperbolóide Tipo-Tempo:**
- Superfície amarela: vetores com n·n = -1
- Folha futura (x⁰ > 0)
- Folha passada (x⁰ < 0)

**3. Vetores Normais:**
- Verde (lime): Orientação futura (T+)
- Vermelho: Orientação passada (T-)
- Escala ajustada para visibilidade
- Labels identificadores

**4. Vetor de Fechamento:**
- Branco: Σ jᵢnᵢ (deve ser ≈ 0)
- Mostra violação de fechamento

#### **Controles Interativos:**

**Configuração:**
- Spin j₀
- Número de trials
- Padrão de orientação:
  - Uniforme (T+) - todos futuros
  - Misto (T+/T-) - 3 futuros, 1 passado
  - Custom patterns disponíveis

**Opções de Display:**
- ☑ Show Light Cone
- ☑ Show Normal Vectors
- ☑ Show Closure Vector
- ☑ Show Hyperboloid

**Controle de Câmera:**
- Slider de elevação (0° - 90°)
- Slider de azimute (-180° - 180°)
- Rotação interativa em tempo real

#### **Interpretação Visual:**

**Configuração Uniforme (T+):**
```
Todos vetores verdes → mesmo cone → podem fechar
det(G) ≈ 10^-60 ✓
```

**Configuração Mista (T+/T-):**
```
Vetores verdes + vermelho → cones opostos → NÃO podem fechar
det(G) ≈ 10^-15 ✗
```

**A obstrução causal é VISÍVEL geometricamente!**

### 📚 **Tab: Theory**
- Resumo teórico completo
- Agora inclui seção sobre visualização 3D
- Equações principais
- Referências

### ✓ **Tab: Basic Check (Table 1)**
Reproduz a Tabela 1 do Apêndice A.5:
- Configuração uniforme (T+)
- Configuração mista (T+/T-)
- Cálculo de det(G)
- Verificação da obstrução causal

### 📈 **Tab: Scaling Analysis**
Análise de escalonamento com spin j:
- Teste múltiplos valores de j
- Verifica persistência da obstrução
- Gráficos log-linear

### 🎯 **Tab: Sensitivity**
Análise de sensibilidade:
- Múltiplas execuções com inicializações aleatórias
- Estatísticas (min, mean, std)
- Histogramas

## 📊 Exemplo de Uso (3D Visualization)

### Passo 1: Configurar
1. Abra tab "3D Visualization"
2. Configure:
   - j₀ = 100
   - trials = 10
   - Orientation = "Mixed (T+/T-)"

### Passo 2: Computar
1. Clique "Compute Configuration"
2. Aguarde cálculo (veja barra de progresso)

### Passo 3: Explorar
1. Use sliders para rotacionar
2. Ative/desative elementos visuais
3. Observe:
   - Vetores verdes (T+) no cone futuro
   - Vetor vermelho (T-) no cone passado
   - Vetor branco (closure) ≠ 0 → obstruction!

### Passo 4: Comparar
1. Mude para "Uniform (T+)"
2. Recompute
3. Observe:
   - Todos vetores verdes
   - Vetor branco ≈ 0 → closure OK!

## 🎨 Interpretação Visual

### Cone de Luz
O cone de luz de Minkowski divide o espaço-tempo em:
- **Interior do cone**: Vetores tipo-tempo (|v⁰| > |v⃗|)
- **Superfície do cone**: Vetores tipo-luz (|v⁰| = |v⃗|)
- **Exterior do cone**: Vetores tipo-espaço (|v⁰| < |v⃗|)

### Hiperbolóide
A condição n·n = -1 define uma hiperbolóide de duas folhas:
- **Folha superior**: Vetores tipo-tempo futuros
- **Folha inferior**: Vetores tipo-tempo passados

### Vetores Normais
Cada vetor normal nᵢ:
- Vive no hiperbolóide (n·n = -1)
- Tem orientação temporal (futuro/passado)
- Contribui para o fechamento com peso jᵢ

### Obstrução Causal
**Geometricamente:**
- Vetores futuros formam um **cone convexo**
- Soma de vetores futuros = vetor futuro
- Não pode cancelar com vetor passado!

**Numericamente:**
- Closure norm ≈ 10^-60 (uniform) → fechamento OK
- Closure norm ≈ 10^-15 (mixed) → fechamento FALHA

## 🎮 Controles de Câmera

### Elevação (0° - 90°)
- 0°: Vista lateral
- 45°: Vista isométrica
- 90°: Vista de cima

### Azimute (-180° - 180°)
- -90°: Vista de trás
- 0°: Vista frontal
- 90°: Vista lateral direita

### Dicas de Visualização:
- **Default (elev=20°, azim=-60°)**: Boa vista geral
- **Para cone futuro**: elev=30°, azim=-45°
- **Para simetria**: elev=20°, azim=0°
- **Vista de topo**: elev=90°, azim=qualquer

## 🔬 Resultados Esperados (3D Viz)

### Uniform (T+):
```
Vetores: 4 verdes (todos no cone futuro)
Closure: ‖Σ jᵢnᵢ‖ ≈ 10^-8 (numericamente ~0)
det(G): ≈ 10^-60
Interpretação: ✓ Geometria Lorentziana real existe
```

### Mixed (T+/T-):
```
Vetores: 3 verdes + 1 vermelho
Closure: ‖Σ jᵢnᵢ‖ ≈ 1.0 (grande!)
det(G): ≈ 10^-15
Interpretação: ✗ Apenas saddle complexo
```

## 📐 Matemática da Visualização

### Projeção
Plotamos (n¹, n², n⁰) em vez de (n⁰, n¹, n²) para melhor visualização.

### Escala
Vetores são escalados por fator 1.5 para visibilidade.

### Cores
```python
Verde (lime): RGB(0, 255, 0)   → Futuro (T+)
Vermelho:     RGB(255, 0, 0)   → Passado (T-)
Branco:       RGB(255, 255, 255) → Closure
Cyan:         Cone futuro
Magenta:      Cone passado
Amarelo:      Hiperbolóide
```

## 🐛 Troubleshooting

### Visualização não aparece
1. Verifique matplotlib backend
2. Tente: `export QT_QPA_PLATFORM=xcb`
3. Reinstale PyQt5

### Performance lenta
1. Reduza `n_points` no código (linha ~1150)
2. Desative hiperbolóide
3. Use menos trials

### Vetores não visíveis
1. Aumente `scale` no código (linha ~1250)
2. Ajuste limites dos eixos
3. Rode mais trials para melhor convergência

## 📊 Comparação v4.1 → v5.0

| Funcionalidade | v4.0 | v4.1 |
|----------------|------|------|
| Basic Check | ✓ | ✓ |
| Scaling | ✓ | ✓ |
| Sensitivity | ✓ | ✓ |
| 3D Viz | Placeholder | **✓ Completo!** |
| Interactive rotation | ✗ | ✓ |
| Light cone | ✗ | ✓ |
| Hyperboloid | ✗ | ✓ |
| Normal vectors | ✗ | ✓ |
| Closure vector | ✗ | ✓ |
| Custom patterns | ✗ | ✓ |

## 🚀 Próximas Versões

### v5.1 (planejado):
- [ ] Exportação 3D para imagens PNG/PDF
- [ ] Animação de rotação automática
- [ ] Múltiplas configurações lado a lado
- [ ] Trajetórias no espaço de configurações

### v6.0 (futuro):
- [ ] Extensão para 2+ vértices
- [ ] Visualização da espuma completa
- [ ] Análise espectral do Hessiano
- [ ] Integração com Jupyter

## 🎓 Referências para 3D Visualization

A visualização 3D implementa conceitos de:

1. **Minkowski Spacetime** (Geometria diferencial)
   - Cone de luz: {x : η(x,x) = 0}
   - Hiperbolóide: {n : η(n,n) = -1}

2. **Closure Constraint** (Apêndice A, Eq. 61)
   - Σ_{b≠a} j_{ab} n_{ab} = 0

3. **Causal Obstruction** (Proposição A.4)
   - Convexidade do cone tipo-tempo

## 📝 Citação

Se usar este código, cite:

```bibtex
@software{guilherme2025spinfoam,
  author = {Guilherme Junior, Mário Sérgio},
  title = {Spin-Foam Euclidean Saddles Analysis Tool v5.0},
  year = {2025},
  note = {3D visualization of causal obstruction}
}
```

## 📧 Contato

Mário Sérgio Guilherme Junior
mario.sergio.guilherme.junior@gmail.com

## 📜 Licença

GPL-3.0 - Use livremente, mantenha atribuição.

---

**✨ Enjoy the 3D visualization! ✨**

_"Seeing is believing - the causal obstruction is now visually obvious."_
