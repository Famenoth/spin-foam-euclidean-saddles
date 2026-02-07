"""
Complete GUI for Spin-Foam Euclidean Saddles Analysis
Author: Mário Sérgio Guilherme Junior

This application provides a complete visualization and analysis interface
for the numerical verification of causally confined Euclidean saddles
in spin-foam quantum gravity.
"""

import sys
import numpy as np
from scipy.optimize import minimize
from numpy.linalg import det, norm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QSpinBox, 
                             QTabWidget, QTextEdit, QGroupBox, QGridLayout,
                             QProgressBar, QComboBox, QDoubleSpinBox, QSlider,
<<<<<<< HEAD
                             QSplitter, QCheckBox, QMessageBox, QFileDialog)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer, QMutex, QCoreApplication
from PyQt5.QtGui import QFont, QColor, QPalette
import time
import traceback
import logging
from contextlib import contextmanager
import platform
import signal
import threading
import json
import csv
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Global configuration for the application"""
    
    # Optimization parameters
    DEFAULT_TRIALS = 10
    MAX_TRIALS = 50
    MIN_TRIALS = 1
    MAX_OPTIMIZATION_ITERATIONS = 1000
    OPTIMIZATION_TOLERANCE = 1e-8
    OPTIMIZATION_TIMEOUT_SECONDS = 30
    
    # Physical parameters
    SPIN_MIN = 10
    SPIN_MAX = 1000
    SPIN_DEFAULT = 100
    
    # Thresholds
    ALPHA_CAUSAL = 1e-3
    BETA_BOUNDARY = 0.1
    LAMBDA_CAUSAL_DEFAULT = 10.0
    LAMBDA_CAUSAL_MIN = 0.1
    LAMBDA_CAUSAL_MAX = 100.0
    
    # UI parameters
    PROGRESS_UPDATE_DELAY_MS = 10
    THREAD_WAIT_TIMEOUT_MS = 1000
    MAIN_WINDOW_WIDTH = 1600
    MAIN_WINDOW_HEIGHT = 1000
    
    # Visualization
    ELEVATION_MIN = 0
    ELEVATION_MAX = 90
    ELEVATION_DEFAULT = 20
    AZIMUTH_MIN = -180
    AZIMUTH_MAX = 180
    AZIMUTH_DEFAULT = -60
    
    # Numerical stability
    GRAM_DET_FLOOR = 1e-100
    ERROR_VALUE_THRESHOLD = 1e9
    
    # Validation limits
    MAX_J_VALUES_COUNT = 20
    MAX_SENSITIVITY_TRIALS = 100
    
    # Export settings
    EXPORT_DATE_FORMAT = "%Y%m%d_%H%M%S"
    EXPORT_JSON_INDENT = 2
    
    # Logging
    LOG_LEVEL = "INFO"
    LOG_FILE = "spinfoam_gui.log"
    LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logging():
    """Configure logging system"""
    logging.basicConfig(
        level=getattr(logging, Config.LOG_LEVEL),
        format=Config.LOG_FORMAT,
        handlers=[
            logging.FileHandler(Config.LOG_FILE),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info("="*80)
    logger.info("Spin-Foam Analysis GUI v5.0 - Starting")
    logger.info("="*80)
    return logger

logger = setup_logging()

# =============================================================================
# TIMEOUT CONTEXT MANAGER
# =============================================================================

class TimeoutError(Exception):
    """Raised when optimization exceeds time limit"""
    pass

if platform.system() != 'Windows':
    # Unix-like systems: use signal.alarm
    @contextmanager
    def time_limit(seconds):
        """Context manager for timeout on Unix systems"""
        def signal_handler(signum, frame):
            raise TimeoutError(f"Optimization exceeded {seconds}s")
        
        old_handler = signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(seconds)
        try:
            yield
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
else:
    # Windows: simplified version (full implementation would require threading)
    @contextmanager
    def time_limit(seconds):
        """Context manager for timeout on Windows (simplified)"""
        # Note: This is a placeholder. Full Windows timeout requires
        # more complex threading implementation.
        logger.warning("Timeout not fully supported on Windows")
        yield
=======
                             QSplitter, QCheckBox, QMessageBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer, QMutex
from PyQt5.QtGui import QFont, QColor, QPalette
import time
import traceback
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d

# =============================================================================
# LORENTZIAN GEOMETRY CORE
# =============================================================================

def lorentz_product(u, v):
<<<<<<< HEAD
    """
    Lorentzian inner product: η_μν = diag(-1, 1, 1, 1)
    
    Parameters
    ----------
    u, v : array_like, shape (4,)
        4-vectors in Minkowski space
        
    Returns
    -------
    float
        Inner product -u⁰v⁰ + u¹v¹ + u²v² + u³v³
    """
    try:
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        if u.shape != (4,) or v.shape != (4,):
            raise ValueError(f"Expected shape (4,), got u: {u.shape}, v: {v.shape}")
        return -u[0]*v[0] + np.dot(u[1:], v[1:])
    except Exception as e:
        logger.error(f"Error in lorentz_product: {e}")
        return 0.0

def lorentz_norm_squared(u):
    """Lorentzian norm squared: u·u"""
=======
    """Lorentzian inner product: η_μν = diag(-1, 1, 1, 1)"""
    try:
        return -u[0]*v[0] + np.dot(u[1:], v[1:])
    except:
        return 0.0

def lorentz_norm_squared(u):
    """Lorentzian norm squared"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    try:
        return lorentz_product(u, u)
    except:
        return 0.0

def project_timelike(u, future=True):
<<<<<<< HEAD
    """
    Project vector onto unit timelike hyperboloid n·n = -1
    
    Parameters
    ----------
    u : array_like, shape (4,)
        Initial 4-vector (only spatial part u[1:] is used)
    future : bool, default=True
        If True, project to future sheet (n⁰ > 0)
        
    Returns
    -------
    ndarray, shape (4,)
        Projected timelike unit vector
    """
    try:
        u = np.asarray(u, dtype=float)
        spatial = u[1:4]
=======
    """Project vector onto unit timelike hyperboloid"""
    try:
        spatial = u[1:]
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        s2 = np.dot(spatial, spatial)
        n0 = np.sqrt(1.0 + s2)
        if not future:
            n0 *= -1.0
<<<<<<< HEAD
        return np.array([n0, spatial[0], spatial[1], spatial[2]])
    except:
=======
        return np.array([n0, *spatial])
    except:
        # Fallback to simple timelike vector
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        if future:
            return np.array([1.0, 0.0, 0.0, 0.0])
        else:
            return np.array([-1.0, 0.0, 0.0, 0.0])

def closure_constraint(normals, spins):
<<<<<<< HEAD
    """
    Compute closure constraint: Σ j_i n_i
    
    Parameters
    ----------
    normals : list of ndarray
        List of 4 normal vectors
    spins : list of float
        List of 4 spins
        
    Returns
    -------
    ndarray, shape (4,)
        Closure vector (should be ~0 for valid geometry)
    """
=======
    """Closure constraint: Σ j_i n_i = 0"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    try:
        return sum(j * n for j, n in zip(spins, normals))
    except:
        return np.zeros(4)

def gram_determinant(normals):
<<<<<<< HEAD
    """
    Compute det(G) where G_ij = n_i · n_j
    
    Parameters
    ----------
    normals : list of ndarray
        List of 4 normal vectors
        
    Returns
    -------
    float
        Determinant of Gram matrix
    """
=======
    """Compute det(G) where G_ij = n_i · n_j"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    try:
        G = np.array([[lorentz_product(ni, nj) for nj in normals] 
                      for ni in normals])
        return det(G)
    except:
        return 0.0

<<<<<<< HEAD
def validate_normal(n, future=True, tolerance=1e-6):
    """Validate that a normal satisfies timelike constraints"""
    try:
        n = np.asarray(n, dtype=float)
        if n.shape != (4,):
            return False
        norm_sq = lorentz_norm_squared(n)
        if abs(norm_sq + 1.0) > tolerance:
            return False
        if future and n[0] <= 0:
            return False
        if not future and n[0] >= 0:
            return False
        return True
    except:
        return False

=======
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
# =============================================================================
# 4-SIMPLEX MODEL
# =============================================================================

class FourSimplex:
    """4-simplex with uniform spins"""
<<<<<<< HEAD
    
    def __init__(self, j0=Config.SPIN_DEFAULT):
        if not (Config.SPIN_MIN <= j0 <= Config.SPIN_MAX):
            raise ValueError(f"j0={j0} out of range [{Config.SPIN_MIN}, {Config.SPIN_MAX}]")
        self.j0 = int(j0)
        self.spins = [self.j0] * 4
        logger.info(f"Created 4-simplex with j0={self.j0}")
    
    def __repr__(self):
        return f"FourSimplex(j0={self.j0})"
=======
    def __init__(self, j0=100):
        self.spins = [j0] * 4
        self.j0 = j0
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d

# =============================================================================
# OBJECTIVE FUNCTION
# =============================================================================

<<<<<<< HEAD
def objective_function(x, spins, causal_pattern, lambda_causal=Config.LAMBDA_CAUSAL_DEFAULT):
    """
    Optimization objective for 4-simplex normal configuration
    
    Minimizes: F = ||closure||² + ||norm violations||² + λ ||causal violations||²
    
    Parameters
    ----------
    x : array_like, shape (16,)
        Flattened array of 4 normal vectors
    spins : list of float
        Spin values
    causal_pattern : list of bool
        Temporal orientations
    lambda_causal : float
        Weight for causal constraints
        
    Returns
    -------
    float
        Objective function value
=======
def objective_function(x, spins, causal_pattern, lambda_causal=10.0):
    """
    Objective function for optimization:
    F = ||closure||² + ||norm violations||² + λ * ||causal violations||²
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    """
    try:
        normals = []
        for i in range(4):
            n = project_timelike(x[4*i:4*i+4], future=causal_pattern[i])
            normals.append(n)

<<<<<<< HEAD
        closure = closure_constraint(normals, spins)
        closure_violation = norm(closure)**2
        
        norm_violation = sum((lorentz_norm_squared(n) + 1.0)**2 for n in normals)
        
=======
        # Closure violation
        closure = closure_constraint(normals, spins)
        closure_violation = norm(closure)**2

        # Norm violation (should be -1 for timelike)
        norm_violation = sum((lorentz_norm_squared(n) + 1.0)**2 for n in normals)

        # Causal violation
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        causal_violation = 0.0
        for n, future in zip(normals, causal_pattern):
            if future:
                causal_violation += max(0.0, -n[0])**2
            else:
                causal_violation += max(0.0, n[0])**2

        return closure_violation + norm_violation + lambda_causal * causal_violation
    except Exception as e:
<<<<<<< HEAD
        logger.error(f"Error in objective_function: {e}")
        return Config.ERROR_VALUE_THRESHOLD

# =============================================================================
# SOLVER WITH TIMEOUT
# =============================================================================

def solve_configuration(simplex, causal_pattern, trials=Config.DEFAULT_TRIALS, 
                       lambda_causal=Config.LAMBDA_CAUSAL_DEFAULT, 
                       callback=None, return_normals=True,
                       timeout=Config.OPTIMIZATION_TIMEOUT_SECONDS):
    """
    Solve for normals satisfying constraints
    
    Parameters
    ----------
    simplex : FourSimplex
        4-simplex model
    causal_pattern : list of bool
        Temporal orientations
    trials : int
        Number of random initializations
    lambda_causal : float
        Weight for causal constraints
    callback : callable or None
        Progress callback
    return_normals : bool
        If False, returns (F_min, None) to save memory
    timeout : int
        Maximum seconds per trial
        
    Returns
    -------
    best_val : float
        Minimum objective value
    best_normals : list of ndarray or None
        Optimal normals
    """
    spins = simplex.spins
    best_val = np.inf
    best_normals = None
    successful_trials = 0
=======
        print(f"Error in objective_function: {e}")
        return 1e10  # Return large value on error

# =============================================================================
# SOLVER
# =============================================================================

def solve_configuration(simplex, causal_pattern, trials=10, 
                       lambda_causal=10.0, callback=None):
    """Solve for normals satisfying constraints"""
    spins = simplex.spins
    best_val = np.inf
    best_normals = None
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    
    for trial in range(trials):
        try:
            x0 = np.random.randn(16)
<<<<<<< HEAD
            
            with time_limit(timeout):
                res = minimize(
                    objective_function, x0,
                    args=(spins, causal_pattern, lambda_causal),
                    method="BFGS",
                    options={
                        "maxiter": Config.MAX_OPTIMIZATION_ITERATIONS,
                        "gtol": Config.OPTIMIZATION_TOLERANCE
                    }
                )
            
            if res.fun < best_val:
                best_val = res.fun
                if return_normals:
                    best_normals = [
                        project_timelike(res.x[4*i:4*i+4], future=causal_pattern[i])
                        for i in range(4)
                    ]
            
            successful_trials += 1
=======
            res = minimize(
                objective_function,
                x0,
                args=(spins, causal_pattern, lambda_causal),
                method="BFGS",
                options={"maxiter": 1000, "gtol": 1e-8}  # Reduced iterations for stability
            )

            if res.fun < best_val:
                best_val = res.fun
                best_normals = [
                    project_timelike(res.x[4*i:4*i+4], future=causal_pattern[i])
                    for i in range(4)
                ]
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            if callback:
                callback(trial + 1, trials)
                
<<<<<<< HEAD
        except TimeoutError:
            logger.warning(f"Trial {trial} timed out after {timeout}s")
            continue
        except Exception as e:
            logger.error(f"Error in trial {trial}: {e}")
            continue
    
    if best_normals is None and return_normals:
        logger.warning("All trials failed, returning fallback normals")
        best_normals = [project_timelike(np.random.randn(4), future=p) 
                       for p in causal_pattern]
        best_val = Config.ERROR_VALUE_THRESHOLD
    
    logger.info(f"Optimization: {successful_trials}/{trials} trials successful, F_min={best_val:.3e}")
    return best_val, best_normals

# =============================================================================
# INPUT VALIDATION
# =============================================================================

def validate_j_values(j_values_str, parent_widget=None):
    """
    Validate and parse j values from user input
    
    Returns
    -------
    list of int or None
        Parsed j values, or None if validation fails
    """
    try:
        j_values = [int(x.strip()) for x in j_values_str.split(',') if x.strip()]
        
        if not j_values:
            if parent_widget:
                QMessageBox.critical(parent_widget, "Invalid Input",
                    "No j values provided.\n\nUse format: 50, 100, 200")
            return None
        
        invalid = [j for j in j_values if not (Config.SPIN_MIN <= j <= Config.SPIN_MAX)]
        if invalid:
            if parent_widget:
                QMessageBox.warning(parent_widget, "Invalid Range",
                    f"Out of range [{Config.SPIN_MIN}, {Config.SPIN_MAX}]:\n{', '.join(map(str, invalid))}")
            return None
        
        if len(j_values) > Config.MAX_J_VALUES_COUNT:
            if parent_widget:
                reply = QMessageBox.question(parent_widget, "Large Computation",
                    f"{len(j_values)} values (max recommended: {Config.MAX_J_VALUES_COUNT}).\n\n" +
                    "This may take a long time. Continue?",
                    QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.No:
                    return None
        
        logger.info(f"Validated j values: {j_values}")
        return j_values
    
    except ValueError as e:
        if parent_widget:
            QMessageBox.critical(parent_widget, "Parse Error",
                f"Invalid format: {e}\n\nUse: 50, 100, 200")
        return None

def validate_spin_parameter(j0, parent_widget=None):
    """Validate single spin parameter"""
    if not (Config.SPIN_MIN <= j0 <= Config.SPIN_MAX):
        if parent_widget:
            QMessageBox.warning(parent_widget, "Invalid Parameter",
                f"j₀={j0} out of range [{Config.SPIN_MIN}, {Config.SPIN_MAX}]")
        return False
    return True

def validate_trials_parameter(trials, parent_widget=None):
    """Validate number of trials"""
    if not (Config.MIN_TRIALS <= trials <= Config.MAX_TRIALS):
        if parent_widget:
            QMessageBox.warning(parent_widget, "Invalid Parameter",
                f"Trials ({trials}) out of range [{Config.MIN_TRIALS}, {Config.MAX_TRIALS}]")
        return False
    
    if trials > 20 and parent_widget:
        reply = QMessageBox.question(parent_widget, "Confirm",
            f"{trials} trials may take several minutes. Continue?",
            QMessageBox.Yes | QMessageBox.No)
        return reply == QMessageBox.Yes
    
    return True

# =============================================================================
# EXPORT UTILITIES
# =============================================================================

def export_results_json(data, filename):
    """Export results to JSON"""
    try:
        with open(filename, 'w') as f:
            json.dump(data, f, indent=Config.EXPORT_JSON_INDENT, default=str)
        logger.info(f"Exported to JSON: {filename}")
        return True
    except Exception as e:
        logger.error(f"JSON export failed: {e}")
        return False

def export_scaling_csv(results, filename):
    """Export scaling analysis to CSV"""
    try:
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['j', 'F_uniform', 'det_G_uniform', 'F_mixed', 'det_G_mixed'])
            for u, m in zip(results['uniform'], results['mixed']):
                writer.writerow([u['j'], u['F'], u['det_G'], m['F'], m['det_G']])
        logger.info(f"Exported to CSV: {filename}")
        return True
    except Exception as e:
        logger.error(f"CSV export failed: {e}")
        return False

# =============================================================================
# COMPUTATION THREAD
=======
        except Exception as e:
            print(f"Error in trial {trial}: {e}")
            continue  # Continue to next trial on error
    
    # Ensure we always return valid normals
    if best_normals is None:
        best_normals = [project_timelike(np.random.randn(4), future=p) 
                       for p in causal_pattern]
        best_val = 1e10
    
    return best_val, best_normals

# =============================================================================
# COMPUTATION THREAD WITH IMPROVED ERROR HANDLING
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
# =============================================================================

class ComputationThread(QThread):
    """Thread for running computations without blocking GUI"""
<<<<<<< HEAD
    progress = pyqtSignal(int, int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
=======
    progress = pyqtSignal(int, int)  # current, total
    finished = pyqtSignal(dict)  # results
    error = pyqtSignal(str)  # error message
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    
    def __init__(self, computation_type, params):
        super().__init__()
        self.computation_type = computation_type
        self.params = params
        self._is_running = True
        self.mutex = QMutex()
        
    def stop(self):
<<<<<<< HEAD
=======
        """Stop the computation"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.mutex.lock()
        self._is_running = False
        self.mutex.unlock()
    
    def is_running(self):
<<<<<<< HEAD
=======
        """Check if thread should continue"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.mutex.lock()
        running = self._is_running
        self.mutex.unlock()
        return running
        
    def run(self):
        try:
            if self.computation_type == "basic":
                self.run_basic_check()
            elif self.computation_type == "scaling":
                self.run_scaling_analysis()
            elif self.computation_type == "sensitivity":
                self.run_sensitivity_analysis()
            elif self.computation_type == "visualization":
                self.run_visualization()
        except Exception as e:
            error_msg = f"Error in {self.computation_type}: {str(e)}\n{traceback.format_exc()}"
<<<<<<< HEAD
            logger.error(error_msg)
=======
            print(error_msg)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            self.error.emit(error_msg)
    
    def run_basic_check(self):
        try:
<<<<<<< HEAD
            j0 = self.params.get('j0', Config.SPIN_DEFAULT)
            trials = self.params.get('trials', Config.DEFAULT_TRIALS)
            lambda_c = self.params.get('lambda_causal', Config.LAMBDA_CAUSAL_DEFAULT)
=======
            j0 = self.params.get('j0', 100)
            trials = self.params.get('trials', 10)
            lambda_c = self.params.get('lambda_causal', 10.0)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            simplex = FourSimplex(j0=j0)
            
            if not self.is_running():
                return
            
<<<<<<< HEAD
            self.progress.emit(0, 2)
            val_u, n_u = solve_configuration(simplex, [True]*4, trials=trials, lambda_causal=lambda_c)
=======
            # Uniform orientation
            self.progress.emit(0, 2)
            val_u, n_u = solve_configuration(
                simplex, [True]*4, trials=trials, lambda_causal=lambda_c
            )
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            if not self.is_running():
                return
            
<<<<<<< HEAD
            self.progress.emit(1, 2)
            val_m, n_m = solve_configuration(simplex, [True, True, True, False], trials=trials, lambda_causal=lambda_c)
=======
            # Mixed orientation
            self.progress.emit(1, 2)
            val_m, n_m = solve_configuration(
                simplex, [True, True, True, False], trials=trials, lambda_causal=lambda_c
            )
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            if not self.is_running():
                return
            
            results = {
                'uniform': {'F': val_u, 'det_G': gram_determinant(n_u), 'normals': n_u},
                'mixed': {'F': val_m, 'det_G': gram_determinant(n_m), 'normals': n_m}
            }
            
            self.progress.emit(2, 2)
            self.finished.emit(results)
            
        except Exception as e:
            self.error.emit(f"Basic check error: {str(e)}")
    
    def run_scaling_analysis(self):
        try:
            j_values = self.params.get('j_values', [50, 100, 200, 500])
            trials = self.params.get('trials', 5)
<<<<<<< HEAD
            lambda_c = self.params.get('lambda_causal', Config.LAMBDA_CAUSAL_DEFAULT)
            
            results = {'j_values': j_values, 'uniform': [], 'mixed': []}
=======
            lambda_c = self.params.get('lambda_causal', 10.0)
            
            results = {'j_values': j_values, 'uniform': [], 'mixed': []}
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            total_steps = len(j_values) * 2
            current_step = 0
            
            for j in j_values:
                if not self.is_running():
                    return
                
                simplex = FourSimplex(j0=j)
                
<<<<<<< HEAD
                self.progress.emit(current_step, total_steps)
                QCoreApplication.processEvents()
                
                try:
                    vu, nu = solve_configuration(simplex, [True]*4, trials=trials, 
                                                lambda_causal=lambda_c, return_normals=True)
                    results['uniform'].append({'j': j, 'F': vu, 'det_G': gram_determinant(nu)})
                except Exception as e:
                    logger.error(f"Error in uniform j={j}: {e}")
                    results['uniform'].append({'j': j, 'F': Config.ERROR_VALUE_THRESHOLD, 'det_G': 0.0})
=======
                # Uniform
                self.progress.emit(current_step, total_steps)
                QThread.msleep(10)  # Small delay to update UI
                
                try:
                    vu, nu = solve_configuration(simplex, [True]*4, trials=trials, 
                                                lambda_causal=lambda_c)
                    results['uniform'].append({
                        'j': j, 'F': vu, 'det_G': gram_determinant(nu)
                    })
                except Exception as e:
                    print(f"Error in uniform j={j}: {e}")
                    results['uniform'].append({
                        'j': j, 'F': 1e10, 'det_G': 0.0
                    })
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
                
                current_step += 1
                
                if not self.is_running():
                    return
                
<<<<<<< HEAD
                self.progress.emit(current_step, total_steps)
                QCoreApplication.processEvents()
                
                try:
                    vm, nm = solve_configuration(simplex, [True, True, True, False], 
                                                trials=trials, lambda_causal=lambda_c, return_normals=True)
                    results['mixed'].append({'j': j, 'F': vm, 'det_G': gram_determinant(nm)})
                except Exception as e:
                    logger.error(f"Error in mixed j={j}: {e}")
                    results['mixed'].append({'j': j, 'F': Config.ERROR_VALUE_THRESHOLD, 'det_G': 0.0})
=======
                # Mixed
                self.progress.emit(current_step, total_steps)
                QThread.msleep(10)
                
                try:
                    vm, nm = solve_configuration(simplex, [True, True, True, False], 
                                                trials=trials, lambda_causal=lambda_c)
                    results['mixed'].append({
                        'j': j, 'F': vm, 'det_G': gram_determinant(nm)
                    })
                except Exception as e:
                    print(f"Error in mixed j={j}: {e}")
                    results['mixed'].append({
                        'j': j, 'F': 1e10, 'det_G': 0.0
                    })
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
                
                current_step += 1
            
            if self.is_running():
<<<<<<< HEAD
                self.progress.emit(total_steps, total_steps)
                QThread.msleep(50)
=======
                # Emit 100% progress before finishing
                self.progress.emit(total_steps, total_steps)
                QThread.msleep(50)  # Brief pause to ensure UI updates
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
                self.finished.emit(results)
                
        except Exception as e:
            self.error.emit(f"Scaling analysis error: {str(e)}")
    
    def run_sensitivity_analysis(self):
        try:
<<<<<<< HEAD
            j0 = self.params.get('j0', Config.SPIN_DEFAULT)
            n_trials = self.params.get('n_trials', 20)
            lambda_c = self.params.get('lambda_causal', Config.LAMBDA_CAUSAL_DEFAULT)
            
            simplex = FourSimplex(j0=j0)
=======
            j0 = self.params.get('j0', 100)
            n_trials = self.params.get('n_trials', 20)
            lambda_c = self.params.get('lambda_causal', 10.0)
            
            simplex = FourSimplex(j0=j0)
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            vals_u, vals_m = [], []
            total_steps = n_trials * 2
            current_step = 0
            
            for i in range(n_trials):
                if not self.is_running():
                    return
                
<<<<<<< HEAD
                self.progress.emit(current_step, total_steps)
                QCoreApplication.processEvents()
                
                try:
                    vu, _ = solve_configuration(simplex, [True]*4, trials=1, 
                                               lambda_causal=lambda_c, return_normals=False)
                    vals_u.append(vu)
                except Exception as e:
                    logger.error(f"Error in uniform trial {i}: {e}")
                    vals_u.append(Config.ERROR_VALUE_THRESHOLD)
=======
                # Uniform
                self.progress.emit(current_step, total_steps)
                QThread.msleep(10)
                
                try:
                    vu, _ = solve_configuration(simplex, [True]*4, trials=1, 
                                               lambda_causal=lambda_c)
                    vals_u.append(vu)
                except Exception as e:
                    print(f"Error in uniform trial {i}: {e}")
                    vals_u.append(1e10)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
                
                current_step += 1
                
                if not self.is_running():
                    return
                
<<<<<<< HEAD
                self.progress.emit(current_step, total_steps)
                QCoreApplication.processEvents()
                
                try:
                    vm, _ = solve_configuration(simplex, [True, True, True, False], 
                                               trials=1, lambda_causal=lambda_c, return_normals=False)
                    vals_m.append(vm)
                except Exception as e:
                    logger.error(f"Error in mixed trial {i}: {e}")
                    vals_m.append(Config.ERROR_VALUE_THRESHOLD)
                
                current_step += 1
            
            results = {'uniform': vals_u, 'mixed': vals_m}
            
            if self.is_running():
                self.progress.emit(total_steps, total_steps)
                QThread.msleep(50)
=======
                # Mixed
                self.progress.emit(current_step, total_steps)
                QThread.msleep(10)
                
                try:
                    vm, _ = solve_configuration(simplex, [True, True, True, False], 
                                               trials=1, lambda_causal=lambda_c)
                    vals_m.append(vm)
                except Exception as e:
                    print(f"Error in mixed trial {i}: {e}")
                    vals_m.append(1e10)
                
                current_step += 1
            
            results = {
                'uniform': vals_u,
                'mixed': vals_m
            }
            
            if self.is_running():
                # Emit 100% progress before finishing
                self.progress.emit(total_steps, total_steps)
                QThread.msleep(50)  # Brief pause to ensure UI updates
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
                self.finished.emit(results)
                
        except Exception as e:
            self.error.emit(f"Sensitivity analysis error: {str(e)}")
    
    def run_visualization(self):
<<<<<<< HEAD
        try:
            j0 = self.params.get('j0', Config.SPIN_DEFAULT)
            trials = self.params.get('trials', Config.DEFAULT_TRIALS)
=======
        """Compute normals for visualization"""
        try:
            j0 = self.params.get('j0', 100)
            trials = self.params.get('trials', 10)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            causal_pattern = self.params.get('causal_pattern', [True]*4)
            
            simplex = FourSimplex(j0=j0)
            
            self.progress.emit(0, 1)
            val, normals = solve_configuration(simplex, causal_pattern, trials=trials)
            
            results = {
                'normals': normals,
                'F': val,
                'det_G': gram_determinant(normals),
                'causal_pattern': causal_pattern,
                'spins': simplex.spins
            }
            
            self.progress.emit(1, 1)
            self.finished.emit(results)
            
        except Exception as e:
            self.error.emit(f"Visualization error: {str(e)}")

# =============================================================================
# MATPLOTLIB CANVAS
# =============================================================================

class MplCanvas(FigureCanvas):
    """Base canvas for matplotlib plots"""
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)

class Mpl3DCanvas(FigureCanvas):
    """Canvas for 3D plots"""
    def __init__(self, parent=None, width=12, height=10, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi, facecolor='#2b2b2b')
        self.axes = self.fig.add_subplot(111, projection='3d', facecolor='#1e1e1e')
<<<<<<< HEAD
        self.fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
=======
        
        # Tight layout for better space usage
        self.fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        super().__init__(self.fig)
        self.setParent(parent)

# =============================================================================
<<<<<<< HEAD
# BASE TAB CLASS (REDUCES CODE DUPLICATION)
# =============================================================================

class BaseAnalysisTab(QWidget):
    """Base class for analysis tabs with common functionality"""
    
    def __init__(self):
        super().__init__()
        self.thread = None
        self.last_results = None
    
    def start_computation(self, thread):
        """Standard protocol for starting computation"""
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.thread.wait()
        
        self.btn_run.setEnabled(False)
        if hasattr(self, 'btn_stop'):
            self.btn_stop.setEnabled(True)
        self.progress.setValue(0)
        if hasattr(self, 'results_text'):
            self.results_text.clear()
        
        self.thread = thread
        self.thread.progress.connect(self.update_progress)
        self.thread.error.connect(self.handle_error)
        self.thread.start()
    
    def finish_computation(self):
        """Standard protocol for finishing computation"""
        self.btn_run.setEnabled(True)
        if hasattr(self, 'btn_stop'):
            self.btn_stop.setEnabled(False)
        if hasattr(self, 'btn_export'):
            self.btn_export.setEnabled(True)
    
    def update_progress(self, current, total):
        """Update progress bar"""
        self.progress.setMaximum(total)
        self.progress.setValue(current)
    
    def handle_error(self, error_msg):
        """Handle computation errors"""
        if hasattr(self, 'results_text'):
            self.results_text.setText(f"ERROR:\n{error_msg}")
        self.finish_computation()
        QMessageBox.warning(self, "Computation Error", 
                          "An error occurred. See results panel for details.")
    
    def export_results(self):
        """Export results (to be overridden by subclasses)"""
        if not self.last_results:
            QMessageBox.information(self, "No Data", "No results to export.")
            return
        
        timestamp = datetime.now().strftime(Config.EXPORT_DATE_FORMAT)
        filename, filter_type = QFileDialog.getSaveFileName(
            self, "Export Results",
            f"results_{timestamp}.json",
            "JSON Files (*.json);;CSV Files (*.csv);;All Files (*)"
        )
        
        if filename:
            if filename.endswith('.json') or filter_type.startswith('JSON'):
                success = export_results_json(self.last_results, filename)
            elif filename.endswith('.csv') or filter_type.startswith('CSV'):
                success = self._export_csv(filename)
            else:
                success = export_results_json(self.last_results, filename)
            
            if success:
                QMessageBox.information(self, "Success", f"Exported to:\n{filename}")
    
    def _export_csv(self, filename):
        """CSV export (to be overridden by subclasses if needed)"""
        QMessageBox.information(self, "Not Implemented", 
                              "CSV export not implemented for this tab.")
        return False

# =============================================================================
# BASIC CHECK TAB
# =============================================================================

class BasicCheckTab(BaseAnalysisTab):
    """Tab for basic verification (Table 1 from paper)"""
    
    def __init__(self):
        super().__init__()
=======
# VISUALIZATION TABS WITH IMPROVED ERROR HANDLING
# =============================================================================

class BasicCheckTab(QWidget):
    """Tab for basic verification (Table 1 from paper)"""
    def __init__(self):
        super().__init__()
        self.thread = None
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Controls
        controls = QGroupBox("Parameters")
        controls_layout = QGridLayout()
        
        controls_layout.addWidget(QLabel("Spin j₀:"), 0, 0)
        self.spin_j0 = QSpinBox()
<<<<<<< HEAD
        self.spin_j0.setRange(Config.SPIN_MIN, Config.SPIN_MAX)
        self.spin_j0.setValue(Config.SPIN_DEFAULT)
=======
        self.spin_j0.setRange(10, 1000)
        self.spin_j0.setValue(100)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls_layout.addWidget(self.spin_j0, 0, 1)
        
        controls_layout.addWidget(QLabel("Trials:"), 1, 0)
        self.spin_trials = QSpinBox()
<<<<<<< HEAD
        self.spin_trials.setRange(Config.MIN_TRIALS, Config.MAX_TRIALS)
        self.spin_trials.setValue(Config.DEFAULT_TRIALS)
=======
        self.spin_trials.setRange(1, 50)
        self.spin_trials.setValue(10)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls_layout.addWidget(self.spin_trials, 1, 1)
        
        controls_layout.addWidget(QLabel("λ (causal):"), 2, 0)
        self.spin_lambda = QDoubleSpinBox()
<<<<<<< HEAD
        self.spin_lambda.setRange(Config.LAMBDA_CAUSAL_MIN, Config.LAMBDA_CAUSAL_MAX)
        self.spin_lambda.setValue(Config.LAMBDA_CAUSAL_DEFAULT)
=======
        self.spin_lambda.setRange(0.1, 100.0)
        self.spin_lambda.setValue(10.0)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls_layout.addWidget(self.spin_lambda, 2, 1)
        
        self.btn_run = QPushButton("Run Analysis")
        self.btn_run.clicked.connect(self.run_analysis)
        controls_layout.addWidget(self.btn_run, 3, 0, 1, 2)
        
<<<<<<< HEAD
        self.btn_export = QPushButton("💾 Export Results")
        self.btn_export.clicked.connect(self.export_results)
        self.btn_export.setEnabled(False)
        controls_layout.addWidget(self.btn_export, 4, 0, 1, 2)
        
=======
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls.setLayout(controls_layout)
        layout.addWidget(controls)
        
        # Progress
        self.progress = QProgressBar()
        layout.addWidget(self.progress)
        
        # Results
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout()
<<<<<<< HEAD
=======
        
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFont(QFont("Courier", 10))
        results_layout.addWidget(self.results_text)
<<<<<<< HEAD
=======
        
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        results_group.setLayout(results_layout)
        layout.addWidget(results_group)
        
        # Visualization
        self.canvas = MplCanvas(self, width=8, height=4)
        layout.addWidget(self.canvas)
        
        self.setLayout(layout)
    
    def run_analysis(self):
<<<<<<< HEAD
        j0 = self.spin_j0.value()
        trials = self.spin_trials.value()
        
        if not validate_spin_parameter(j0, self):
            return
        if not validate_trials_parameter(trials, self):
            return
        
        params = {
            'j0': j0,
            'trials': trials,
            'lambda_causal': self.spin_lambda.value()
        }
        
        thread = ComputationThread("basic", params)
        thread.finished.connect(self.display_results)
        self.start_computation(thread)
    
    def display_results(self, results):
        try:
            self.last_results = results
            
            # Text
            text = "=" * 60 + "\n"
            text += "NUMERICAL CONSISTENCY CHECK (Appendix A.5)\n"
            text += "=" * 60 + "\n\n"
            text += f"Parameters: j₀={self.spin_j0.value()}, trials={self.spin_trials.value()}, λ={self.spin_lambda.value():.1f}\n\n"
            text += f"Uniform (T+):   F={results['uniform']['F']:.3e}, |det(G)|={abs(results['uniform']['det_G']):.3e}\n"
            text += f"Mixed (T+/T-):  F={results['mixed']['F']:.3e}, |det(G)|={abs(results['mixed']['det_G']):.3e}\n\n"
            
            det_ratio = abs(results['mixed']['det_G']) / max(abs(results['uniform']['det_G']), Config.GRAM_DET_FLOOR)
            text += f"det(G) ratio (mixed/uniform): {det_ratio:.2e}\n\n"
            
            if det_ratio > 1e10:
                text += "✓ Strong numerical evidence for causal obstruction\n"
            else:
                text += "⚠ Weak evidence - may need parameter tuning\n"
            
            text += "\n" + "=" * 60
=======
        # Stop previous thread if running
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.thread.wait()
        
        self.btn_run.setEnabled(False)
        self.progress.setValue(0)
        self.results_text.clear()
        
        params = {
            'j0': self.spin_j0.value(),
            'trials': self.spin_trials.value(),
            'lambda_causal': self.spin_lambda.value()
        }
        
        self.thread = ComputationThread("basic", params)
        self.thread.progress.connect(self.update_progress)
        self.thread.finished.connect(self.display_results)
        self.thread.error.connect(self.handle_error)
        self.thread.start()
    
    def update_progress(self, current, total):
        self.progress.setMaximum(total)
        self.progress.setValue(current)
    
    def handle_error(self, error_msg):
        self.results_text.setText(f"ERROR:\n{error_msg}")
        self.btn_run.setEnabled(True)
        QMessageBox.warning(self, "Computation Error", 
                          "An error occurred during computation. See results panel for details.")
    
    def display_results(self, results):
        try:
            # Text results
            text = "=" * 60 + "\n"
            text += "NUMERICAL CONSISTENCY CHECK (Appendix A.5)\n"
            text += "=" * 60 + "\n\n"
            
            text += f"Parameters: j₀ = {self.spin_j0.value()}, "
            text += f"trials = {self.spin_trials.value()}, "
            text += f"λ = {self.spin_lambda.value():.1f}\n\n"
            
            text += "Uniform Orientation (T+):\n"
            text += f"  F_min     = {results['uniform']['F']:.3e}\n"
            text += f"  |det(G)|  = {abs(results['uniform']['det_G']):.3e}\n\n"
            
            text += "Mixed Orientation (T+/T-):\n"
            text += f"  F_min     = {results['mixed']['F']:.3e}\n"
            text += f"  |det(G)|  = {abs(results['mixed']['det_G']):.3e}\n\n"
            
            text += "=" * 60 + "\n"
            text += "INTERPRETATION:\n"
            text += "=" * 60 + "\n"
            
            det_ratio = abs(results['mixed']['det_G']) / max(abs(results['uniform']['det_G']), 1e-100)
            text += f"\ndet(G) ratio (mixed/uniform): {det_ratio:.2e}\n\n"
            
            if det_ratio > 1e10:
                text += "✓ Strong numerical evidence for causal obstruction\n"
                text += "✓ Mixed orientation prevents vectorial closure\n"
                text += "✓ Consistent with Proposition A.4 (analytical argument)\n"
            else:
                text += "⚠ Weak evidence - may need parameter tuning\n"
            
            text += "\n" + "=" * 60 + "\n"
            text += "This is a PROOF-OF-CONCEPT consistency check.\n"
            text += "See paper Section A.5 for complete discussion.\n"
            text += "=" * 60 + "\n"
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            self.results_text.setText(text)
            
            # Plot
            self.canvas.axes.clear()
<<<<<<< HEAD
            categories = ['Uniform\n(T+)', 'Mixed\n(T+/T-)']
            det_values = [max(abs(results['uniform']['det_G']), Config.GRAM_DET_FLOOR), 
                         max(abs(results['mixed']['det_G']), Config.GRAM_DET_FLOOR)]
            bars = self.canvas.axes.bar(categories, det_values, color=['green', 'red'])
            self.canvas.axes.set_ylabel('|det(G)|', fontsize=12, fontweight='bold')
            self.canvas.axes.set_title('Gram Determinant: Causal Obstruction', fontsize=14, fontweight='bold')
            self.canvas.axes.set_yscale('log')
            self.canvas.axes.grid(True, alpha=0.3)
            
            for bar, val in zip(bars, det_values):
                height = bar.get_height()
                self.canvas.axes.text(bar.get_x() + bar.get_width()/2., height,
                                     f'{val:.1e}', ha='center', va='bottom', fontweight='bold')
=======
            
            categories = ['Uniform\n(T+)', 'Mixed\n(T+/T-)']
            det_values = [max(abs(results['uniform']['det_G']), 1e-100), 
                         max(abs(results['mixed']['det_G']), 1e-100)]
            
            bars = self.canvas.axes.bar(categories, det_values, color=['green', 'red'])
            self.canvas.axes.set_ylabel('|det(G)|', fontsize=12, fontweight='bold')
            self.canvas.axes.set_title('Gram Determinant: Causal Obstruction', 
                                       fontsize=14, fontweight='bold')
            self.canvas.axes.set_yscale('log')
            self.canvas.axes.grid(True, alpha=0.3)
            
            # Add value labels
            for bar, val in zip(bars, det_values):
                height = bar.get_height()
                self.canvas.axes.text(bar.get_x() + bar.get_width()/2., height,
                                     f'{val:.1e}',
                                     ha='center', va='bottom', fontweight='bold')
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            self.canvas.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            self.handle_error(f"Error displaying results: {str(e)}")
        finally:
<<<<<<< HEAD
            self.finish_computation()

# =============================================================================
# SCALING TAB
# =============================================================================

class ScalingTab(BaseAnalysisTab):
    """Tab for scaling analysis with spin j"""
    
    def __init__(self):
        super().__init__()
=======
            self.btn_run.setEnabled(True)

class ScalingTab(QWidget):
    """Tab for scaling analysis with spin j"""
    def __init__(self):
        super().__init__()
        self.thread = None
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Controls
        controls = QGroupBox("Parameters")
        controls_layout = QGridLayout()
        
        controls_layout.addWidget(QLabel("j values (comma-separated):"), 0, 0)
        self.j_values_input = QTextEdit()
        self.j_values_input.setMaximumHeight(60)
        self.j_values_input.setText("50, 100, 200, 500")
        controls_layout.addWidget(self.j_values_input, 0, 1)
        
        controls_layout.addWidget(QLabel("Trials per j:"), 1, 0)
        self.spin_trials = QSpinBox()
<<<<<<< HEAD
        self.spin_trials.setRange(Config.MIN_TRIALS, 20)
=======
        self.spin_trials.setRange(1, 20)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.spin_trials.setValue(5)
        controls_layout.addWidget(self.spin_trials, 1, 1)
        
        self.btn_run = QPushButton("Run Scaling Analysis")
        self.btn_run.clicked.connect(self.run_analysis)
        controls_layout.addWidget(self.btn_run, 2, 0, 1, 2)
        
        self.btn_stop = QPushButton("Stop Analysis")
        self.btn_stop.clicked.connect(self.stop_analysis)
        self.btn_stop.setEnabled(False)
        self.btn_stop.setStyleSheet("background-color: #a04040;")
        controls_layout.addWidget(self.btn_stop, 3, 0, 1, 2)
        
<<<<<<< HEAD
        self.btn_export = QPushButton("💾 Export Results")
        self.btn_export.clicked.connect(self.export_results)
        self.btn_export.setEnabled(False)
        controls_layout.addWidget(self.btn_export, 4, 0, 1, 2)
        
=======
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls.setLayout(controls_layout)
        layout.addWidget(controls)
        
        # Progress
        self.progress = QProgressBar()
        layout.addWidget(self.progress)
        
        # Results
        splitter = QSplitter(Qt.Horizontal)
<<<<<<< HEAD
=======
        
        # Table
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFont(QFont("Courier", 9))
        splitter.addWidget(self.results_text)
        
<<<<<<< HEAD
=======
        # Plot
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.canvas = MplCanvas(self, width=8, height=6)
        splitter.addWidget(self.canvas)
        
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)
<<<<<<< HEAD
=======
        
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        layout.addWidget(splitter)
        
        self.setLayout(layout)
    
    def run_analysis(self):
<<<<<<< HEAD
        j_values = validate_j_values(self.j_values_input.toPlainText(), self)
        if not j_values:
            return
        
        trials = self.spin_trials.value()
        if not validate_trials_parameter(trials, self):
            return
        
        params = {'j_values': j_values, 'trials': trials, 'lambda_causal': Config.LAMBDA_CAUSAL_DEFAULT}
        thread = ComputationThread("scaling", params)
        thread.finished.connect(self.display_results)
        self.start_computation(thread)
=======
        try:
            j_values = [int(x.strip()) for x in self.j_values_input.toPlainText().split(',')]
            if not j_values:
                raise ValueError("No j values provided")
        except Exception as e:
            self.results_text.setText(f"Error: Invalid j values format\n{str(e)}")
            return
        
        # Stop previous thread if running
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.thread.wait()
        
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.progress.setValue(0)
        self.results_text.clear()
        
        params = {
            'j_values': j_values,
            'trials': self.spin_trials.value(),
            'lambda_causal': 10.0
        }
        
        self.thread = ComputationThread("scaling", params)
        self.thread.progress.connect(self.update_progress)
        self.thread.finished.connect(self.display_results)
        self.thread.error.connect(self.handle_error)
        self.thread.start()
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    
    def stop_analysis(self):
        if self.thread and self.thread.isRunning():
            self.thread.stop()
<<<<<<< HEAD
            self.results_text.append("\n\n*** STOPPED BY USER ***\n")
            self.finish_computation()
    
    def display_results(self, results):
        try:
            self.last_results = results
            
=======
            self.results_text.append("\n\n*** ANALYSIS STOPPED BY USER ***\n")
            self.btn_stop.setEnabled(False)
            self.btn_run.setEnabled(True)
    
    def update_progress(self, current, total):
        self.progress.setMaximum(total)
        self.progress.setValue(current)
    
    def handle_error(self, error_msg):
        self.results_text.setText(f"ERROR:\n{error_msg}")
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        QMessageBox.warning(self, "Computation Error", 
                          "An error occurred during computation. See results panel for details.")
    
    def display_results(self, results):
        try:
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            # Table
            text = "=" * 80 + "\n"
            text += "SCALING ANALYSIS WITH SPIN j\n"
            text += "=" * 80 + "\n\n"
<<<<<<< HEAD
            text += f"{'j':>6} | {'F_uniform':>12} | {'det(G)_uniform':>15} | {'F_mixed':>12} | {'det(G)_mixed':>15}\n"
            text += "-" * 80 + "\n"
            
            for u, m in zip(results['uniform'], results['mixed']):
                text += f"{u['j']:6d} | {u['F']:12.2e} | {abs(u['det_G']):15.2e} | {m['F']:12.2e} | {abs(m['det_G']):15.2e}\n"
            
            text += "\n" + "=" * 80 + "\n"
            text += "The obstruction persists across different spin scales.\n"
            text += "=" * 80
=======
            
            text += f"{'j':>6} | {'F_uniform':>12} | {'det(G)_uniform':>15} | "
            text += f"{'F_mixed':>12} | {'det(G)_mixed':>15}\n"
            text += "-" * 80 + "\n"
            
            for u, m in zip(results['uniform'], results['mixed']):
                text += f"{u['j']:6d} | {u['F']:12.2e} | {abs(u['det_G']):15.2e} | "
                text += f"{m['F']:12.2e} | {abs(m['det_G']):15.2e}\n"
            
            text += "\n" + "=" * 80 + "\n"
            text += "OBSERVATION:\n"
            text += "The obstruction persists across different spin scales,\n"
            text += "supporting the genericity claim in Proposition A.4.\n"
            text += "=" * 80 + "\n"
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            self.results_text.setText(text)
            
            # Plot
            self.canvas.axes.clear()
<<<<<<< HEAD
            j_vals = results['j_values']
            det_u = [max(abs(u['det_G']), Config.GRAM_DET_FLOOR) for u in results['uniform']]
            det_m = [max(abs(m['det_G']), Config.GRAM_DET_FLOOR) for m in results['mixed']]
            
            self.canvas.axes.semilogy(j_vals, det_u, 'o-', label='Uniform (T+)', linewidth=2, markersize=8, color='green')
            self.canvas.axes.semilogy(j_vals, det_m, 's-', label='Mixed (T+/T-)', linewidth=2, markersize=8, color='red')
            
            self.canvas.axes.set_xlabel('Spin j₀', fontsize=12, fontweight='bold')
            self.canvas.axes.set_ylabel('|det(G)|', fontsize=12, fontweight='bold')
            self.canvas.axes.set_title('Scaling of Gram Determinant', fontsize=14, fontweight='bold')
            self.canvas.axes.legend(fontsize=11)
            self.canvas.axes.grid(True, alpha=0.3)
=======
            
            j_vals = results['j_values']
            det_u = [max(abs(u['det_G']), 1e-100) for u in results['uniform']]
            det_m = [max(abs(m['det_G']), 1e-100) for m in results['mixed']]
            
            self.canvas.axes.semilogy(j_vals, det_u, 'o-', label='Uniform (T+)', 
                                     linewidth=2, markersize=8, color='green')
            self.canvas.axes.semilogy(j_vals, det_m, 's-', label='Mixed (T+/T-)', 
                                     linewidth=2, markersize=8, color='red')
            
            self.canvas.axes.set_xlabel('Spin j₀', fontsize=12, fontweight='bold')
            self.canvas.axes.set_ylabel('|det(G)|', fontsize=12, fontweight='bold')
            self.canvas.axes.set_title('Scaling of Gram Determinant with Spin', 
                                       fontsize=14, fontweight='bold')
            self.canvas.axes.legend(fontsize=11)
            self.canvas.axes.grid(True, alpha=0.3)
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            self.canvas.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            self.handle_error(f"Error displaying results: {str(e)}")
        finally:
<<<<<<< HEAD
            self.finish_computation()
    
    def _export_csv(self, filename):
        """Export scaling results to CSV"""
        return export_scaling_csv(self.last_results, filename)

# =============================================================================
# SENSITIVITY TAB
# =============================================================================

class SensitivityTab(BaseAnalysisTab):
    """Tab for sensitivity analysis (multiple trials)"""
    
    def __init__(self):
        super().__init__()
=======
            self.btn_run.setEnabled(True)
            self.btn_stop.setEnabled(False)

class SensitivityTab(QWidget):
    """Tab for sensitivity analysis (multiple trials)"""
    def __init__(self):
        super().__init__()
        self.thread = None
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Controls
        controls = QGroupBox("Parameters")
        controls_layout = QGridLayout()
        
        controls_layout.addWidget(QLabel("Spin j₀:"), 0, 0)
        self.spin_j0 = QSpinBox()
<<<<<<< HEAD
        self.spin_j0.setRange(Config.SPIN_MIN, Config.SPIN_MAX)
        self.spin_j0.setValue(Config.SPIN_DEFAULT)
=======
        self.spin_j0.setRange(10, 1000)
        self.spin_j0.setValue(100)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls_layout.addWidget(self.spin_j0, 0, 1)
        
        controls_layout.addWidget(QLabel("Number of trials:"), 1, 0)
        self.spin_trials = QSpinBox()
<<<<<<< HEAD
        self.spin_trials.setRange(5, Config.MAX_SENSITIVITY_TRIALS)
=======
        self.spin_trials.setRange(5, 100)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.spin_trials.setValue(20)
        controls_layout.addWidget(self.spin_trials, 1, 1)
        
        self.btn_run = QPushButton("Run Sensitivity Analysis")
        self.btn_run.clicked.connect(self.run_analysis)
        controls_layout.addWidget(self.btn_run, 2, 0, 1, 2)
        
        self.btn_stop = QPushButton("Stop Analysis")
        self.btn_stop.clicked.connect(self.stop_analysis)
        self.btn_stop.setEnabled(False)
        self.btn_stop.setStyleSheet("background-color: #a04040;")
        controls_layout.addWidget(self.btn_stop, 3, 0, 1, 2)
        
<<<<<<< HEAD
        self.btn_export = QPushButton("💾 Export Results")
        self.btn_export.clicked.connect(self.export_results)
        self.btn_export.setEnabled(False)
        controls_layout.addWidget(self.btn_export, 4, 0, 1, 2)
        
=======
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls.setLayout(controls_layout)
        layout.addWidget(controls)
        
        # Progress
        self.progress = QProgressBar()
        layout.addWidget(self.progress)
        
        # Results
        splitter = QSplitter(Qt.Horizontal)
        
        # Stats
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFont(QFont("Courier", 10))
        splitter.addWidget(self.results_text)
        
        # Plot
        self.canvas = MplCanvas(self, width=8, height=6)
        splitter.addWidget(self.canvas)
        
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)
        
        layout.addWidget(splitter)
        
        self.setLayout(layout)
    
    def run_analysis(self):
<<<<<<< HEAD
        j0 = self.spin_j0.value()
        trials = self.spin_trials.value()
        
        if not validate_spin_parameter(j0, self):
            return
        if not validate_trials_parameter(trials, self):
            return
        
        params = {
            'j0': j0,
            'n_trials': trials,
            'lambda_causal': Config.LAMBDA_CAUSAL_DEFAULT
        }
        
        thread = ComputationThread("sensitivity", params)
        thread.finished.connect(self.display_results)
        self.start_computation(thread)
=======
        # Stop previous thread if running
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.thread.wait()
        
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.progress.setValue(0)
        self.results_text.clear()
        
        params = {
            'j0': self.spin_j0.value(),
            'n_trials': self.spin_trials.value(),
            'lambda_causal': 10.0
        }
        
        self.thread = ComputationThread("sensitivity", params)
        self.thread.progress.connect(self.update_progress)
        self.thread.finished.connect(self.display_results)
        self.thread.error.connect(self.handle_error)
        self.thread.start()
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    
    def stop_analysis(self):
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.results_text.append("\n\n*** ANALYSIS STOPPED BY USER ***\n")
<<<<<<< HEAD
            self.finish_computation()
    
    def display_results(self, results):
        try:
            self.last_results = results
            
            # Filter out error values
            vals_u = [v for v in results['uniform'] if v < Config.ERROR_VALUE_THRESHOLD]
            vals_m = [v for v in results['mixed'] if v < Config.ERROR_VALUE_THRESHOLD]
            
            if not vals_u or not vals_m:
                self.results_text.setText("ERROR: All trials failed. Try reducing trials or adjusting parameters.")
                self.finish_computation()
                return
            
            # Convert to numpy arrays for statistics
            vals_u_array = np.array(vals_u)
            vals_m_array = np.array(vals_m)
            
=======
            self.btn_stop.setEnabled(False)
            self.btn_run.setEnabled(True)
    
    def update_progress(self, current, total):
        self.progress.setMaximum(total)
        self.progress.setValue(current)
    
    def handle_error(self, error_msg):
        self.results_text.setText(f"ERROR:\n{error_msg}")
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        QMessageBox.warning(self, "Computation Error", 
                          "An error occurred during computation. See results panel for details.")
    
    def display_results(self, results):
        try:
            vals_u = results['uniform']
            vals_m = results['mixed']
            
            # Filter out error values
            vals_u = [v for v in vals_u if v < 1e9]
            vals_m = [v for v in vals_m if v < 1e9]
            
            if not vals_u or not vals_m:
                self.results_text.setText("ERROR: All trials failed. Try reducing number of trials or adjusting parameters.")
                return
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            # Statistics
            text = "=" * 60 + "\n"
            text += "SENSITIVITY ANALYSIS\n"
            text += "=" * 60 + "\n\n"
<<<<<<< HEAD
=======
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            text += f"Configuration: j₀ = {self.spin_j0.value()}\n"
            text += f"Successful trials: {len(vals_u)} uniform, {len(vals_m)} mixed\n\n"
            
            text += "Uniform Orientation (T+):\n"
<<<<<<< HEAD
            text += f"  min(F)  = {np.min(vals_u_array):.3e}\n"
            text += f"  mean(F) = {np.mean(vals_u_array):.3e}\n"
            text += f"  std(F)  = {np.std(vals_u_array):.3e}\n\n"
            
            text += "Mixed Orientation (T+/T-):\n"
            text += f"  min(F)  = {np.min(vals_m_array):.3e}\n"
            text += f"  mean(F) = {np.mean(vals_m_array):.3e}\n"
            text += f"  std(F)  = {np.std(vals_m_array):.3e}\n\n"
=======
            text += f"  min(F)  = {min(vals_u):.3e}\n"
            text += f"  mean(F) = {np.mean(vals_u):.3e}\n"
            text += f"  std(F)  = {np.std(vals_u):.3e}\n\n"
            
            text += "Mixed Orientation (T+/T-):\n"
            text += f"  min(F)  = {min(vals_m):.3e}\n"
            text += f"  mean(F) = {np.mean(vals_m):.3e}\n"
            text += f"  std(F)  = {np.std(vals_m):.3e}\n\n"
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            text += "=" * 60 + "\n"
            text += "INTERPRETATION:\n"
            text += "=" * 60 + "\n"
            
<<<<<<< HEAD
            if np.min(vals_m_array) > 10 * np.max(vals_u_array):
                text += "✓ Mixed configuration consistently fails to close\n"
                text += "✓ Uniform configuration achieves closure reliably\n"
                text += "✓ Results stable across random initializations\n"
=======
            if min(vals_m) > 10 * max(vals_u):
                text += "✓ Mixed configuration consistently fails to close\n"
                text += "✓ Uniform configuration achieves closure reliably\n"
                text += "✓ Results are stable across multiple random initializations\n"
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            else:
                text += "⚠ Results show overlap - may indicate:\n"
                text += "  - Need for more trials\n"
                text += "  - Parameter sensitivity\n"
                text += "  - Numerical conditioning issues\n"
            
<<<<<<< HEAD
            text += "\n" + "=" * 60
=======
            text += "\n" + "=" * 60 + "\n"
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            self.results_text.setText(text)
            
            # Histogram plot
            self.canvas.axes.clear()
            
            self.canvas.axes.hist(vals_u, bins=15, alpha=0.7, label='Uniform (T+)', 
                                 color='green', edgecolor='black')
            self.canvas.axes.hist(vals_m, bins=15, alpha=0.7, label='Mixed (T+/T-)', 
                                 color='red', edgecolor='black')
            
            self.canvas.axes.set_xlabel('Objective Function Value F', fontsize=12, 
                                       fontweight='bold')
            self.canvas.axes.set_ylabel('Frequency', fontsize=12, fontweight='bold')
            self.canvas.axes.set_title('Distribution of Optimization Results', 
                                       fontsize=14, fontweight='bold')
            self.canvas.axes.legend(fontsize=11)
            self.canvas.axes.grid(True, alpha=0.3, axis='y')
            
            # Add vertical lines for means
<<<<<<< HEAD
            mean_u = np.mean(vals_u_array)
            mean_m = np.mean(vals_m_array)
            
            self.canvas.axes.axvline(mean_u, color='darkgreen', 
                                    linestyle='--', linewidth=2, label=f'Mean (Uniform): {mean_u:.2e}')
            self.canvas.axes.axvline(mean_m, color='darkred', 
                                    linestyle='--', linewidth=2, label=f'Mean (Mixed): {mean_m:.2e}')
=======
            self.canvas.axes.axvline(np.mean(vals_u), color='darkgreen', 
                                    linestyle='--', linewidth=2, label='Mean (Uniform)')
            self.canvas.axes.axvline(np.mean(vals_m), color='darkred', 
                                    linestyle='--', linewidth=2, label='Mean (Mixed)')
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            
            self.canvas.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
<<<<<<< HEAD
            error_msg = f"Error displaying results: {str(e)}\n{traceback.format_exc()}"
            logger.error(error_msg)
            self.handle_error(error_msg)
        finally:
            self.finish_computation()

# =============================================================================
# VISUALIZATION TAB (continued in next message due to length)
# =============================================================================

class VisualizationTab(QWidget):
    """Tab for 3D visualization"""
    
=======
            self.handle_error(f"Error displaying results: {str(e)}")
        finally:
            self.btn_run.setEnabled(True)
            self.btn_stop.setEnabled(False)

# [CONTINUAÇÃO NO PRÓXIMO ARQUIVO - código muito longo]
# Esta é a parte 2 - adicione ao final do arquivo anterior

class VisualizationTab(QWidget):
    """Tab for 3D visualization of normals and light cone"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    def __init__(self):
        super().__init__()
        self.current_results = None
        self.fullscreen_mode = False
<<<<<<< HEAD
=======
        self.left_panel_ref = None
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.thread = None
        self.init_ui()
        
    def init_ui(self):
<<<<<<< HEAD
        main_layout = QVBoxLayout()
=======
        # Main layout with splitter for better space management
        main_layout = QVBoxLayout()
        
        # Header
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        header = QLabel("3D Minkowski Light Cone Visualization")
        header.setFont(QFont("Arial", 14, QFont.Bold))
        header.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(header)
        
<<<<<<< HEAD
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        
=======
        # Create splitter for controls and visualization
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel: Controls
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        
        # Configuration controls
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls = QGroupBox("Configuration")
        controls_layout = QGridLayout()
        
        controls_layout.addWidget(QLabel("Spin j₀:"), 0, 0)
        self.spin_j0 = QSpinBox()
<<<<<<< HEAD
        self.spin_j0.setRange(Config.SPIN_MIN, Config.SPIN_MAX)
        self.spin_j0.setValue(Config.SPIN_DEFAULT)
=======
        self.spin_j0.setRange(10, 1000)
        self.spin_j0.setValue(100)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls_layout.addWidget(self.spin_j0, 0, 1)
        
        controls_layout.addWidget(QLabel("Trials:"), 1, 0)
        self.spin_trials = QSpinBox()
        self.spin_trials.setRange(1, 20)
<<<<<<< HEAD
        self.spin_trials.setValue(Config.DEFAULT_TRIALS)
=======
        self.spin_trials.setValue(10)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        controls_layout.addWidget(self.spin_trials, 1, 1)
        
        controls_layout.addWidget(QLabel("Orientation:"), 2, 0)
        self.combo_orientation = QComboBox()
        self.combo_orientation.addItems([
<<<<<<< HEAD
            "Uniform (T+)", "Mixed (T+/T-)", "Custom: T+,T+,T+,T-",
            "Custom: T+,T+,T-,T-", "Custom: T+,T-,T+,T-"
=======
            "Uniform (T+)",
            "Mixed (T+/T-)",
            "Custom: T+,T+,T+,T-",
            "Custom: T+,T+,T-,T-",
            "Custom: T+,T-,T+,T-"
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        ])
        controls_layout.addWidget(self.combo_orientation, 2, 1)
        
        self.btn_compute = QPushButton("Compute Configuration")
        self.btn_compute.clicked.connect(self.compute_configuration)
        controls_layout.addWidget(self.btn_compute, 3, 0, 1, 2)
        
        self.btn_fullscreen = QPushButton("🔍 Maximize 3D View")
        self.btn_fullscreen.clicked.connect(self.toggle_fullscreen)
        self.btn_fullscreen.setStyleSheet("background-color: #146C94;")
        controls_layout.addWidget(self.btn_fullscreen, 4, 0, 1, 2)
        
        controls.setLayout(controls_layout)
        left_layout.addWidget(controls)
        
<<<<<<< HEAD
        self.progress = QProgressBar()
        left_layout.addWidget(self.progress)
        
        # Display options
=======
        # Progress
        self.progress = QProgressBar()
        left_layout.addWidget(self.progress)
        
        # Visualization controls
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        viz_controls = QGroupBox("Display Options")
        viz_layout = QGridLayout()
        
        self.check_cone = QCheckBox("Show Light Cone")
        self.check_cone.setChecked(True)
        self.check_cone.stateChanged.connect(self.update_plot)
        viz_layout.addWidget(self.check_cone, 0, 0, 1, 2)
        
        self.check_normals = QCheckBox("Show Normal Vectors")
        self.check_normals.setChecked(True)
        self.check_normals.stateChanged.connect(self.update_plot)
        viz_layout.addWidget(self.check_normals, 1, 0, 1, 2)
        
        self.check_closure = QCheckBox("Show Closure Vector")
        self.check_closure.setChecked(True)
        self.check_closure.stateChanged.connect(self.update_plot)
        viz_layout.addWidget(self.check_closure, 2, 0, 1, 2)
        
        self.check_hyperboloid = QCheckBox("Show Hyperboloid")
        self.check_hyperboloid.setChecked(True)
        self.check_hyperboloid.stateChanged.connect(self.update_plot)
        viz_layout.addWidget(self.check_hyperboloid, 3, 0, 1, 2)
        
        viz_layout.addWidget(QLabel("Elevation:"), 4, 0)
<<<<<<< HEAD
        self.label_elev = QLabel(f"{Config.ELEVATION_DEFAULT}°")
=======
        self.label_elev = QLabel("20°")
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.label_elev.setAlignment(Qt.AlignRight)
        viz_layout.addWidget(self.label_elev, 4, 1)
        
        self.slider_elev = QSlider(Qt.Horizontal)
<<<<<<< HEAD
        self.slider_elev.setRange(Config.ELEVATION_MIN, Config.ELEVATION_MAX)
        self.slider_elev.setValue(Config.ELEVATION_DEFAULT)
=======
        self.slider_elev.setRange(0, 90)
        self.slider_elev.setValue(20)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.slider_elev.valueChanged.connect(self.update_view)
        viz_layout.addWidget(self.slider_elev, 5, 0, 1, 2)
        
        viz_layout.addWidget(QLabel("Azimuth:"), 6, 0)
<<<<<<< HEAD
        self.label_azim = QLabel(f"{Config.AZIMUTH_DEFAULT}°")
=======
        self.label_azim = QLabel("-60°")
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.label_azim.setAlignment(Qt.AlignRight)
        viz_layout.addWidget(self.label_azim, 6, 1)
        
        self.slider_azim = QSlider(Qt.Horizontal)
<<<<<<< HEAD
        self.slider_azim.setRange(Config.AZIMUTH_MIN, Config.AZIMUTH_MAX)
        self.slider_azim.setValue(Config.AZIMUTH_DEFAULT)
=======
        self.slider_azim.setRange(-180, 180)
        self.slider_azim.setValue(-60)
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.slider_azim.valueChanged.connect(self.update_view)
        viz_layout.addWidget(self.slider_azim, 7, 0, 1, 2)
        
        viz_controls.setLayout(viz_layout)
        left_layout.addWidget(viz_controls)
        
<<<<<<< HEAD
        # Info
=======
        # Info panel
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        info_group = QGroupBox("Configuration Info")
        info_layout = QVBoxLayout()
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setFont(QFont("Courier", 9))
        info_layout.addWidget(self.info_text)
        info_group.setLayout(info_layout)
        left_layout.addWidget(info_group)
        
<<<<<<< HEAD
        left_layout.addStretch()
        left_panel.setLayout(left_layout)
        self.left_panel_ref = left_panel
        splitter.addWidget(left_panel)
        
        # 3D Canvas
        self.canvas = Mpl3DCanvas(self, width=12, height=10)
        splitter.addWidget(self.canvas)
        
        self.splitter_ref = splitter
=======
        # Add stretch to push everything to top
        left_layout.addStretch()
        
        left_panel.setLayout(left_layout)
        self.left_panel_ref = left_panel  # Store reference
        splitter.addWidget(left_panel)
        
        # Right panel: 3D Visualization (MUCH LARGER)
        self.canvas = Mpl3DCanvas(self, width=12, height=10)
        splitter.addWidget(self.canvas)
        
        # Set splitter sizes: 25% controls, 75% visualization
        self.splitter_ref = splitter  # Store reference
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)
        
        main_layout.addWidget(splitter)
<<<<<<< HEAD
        self.setLayout(main_layout)
        self.plot_empty_cone()
    
    def toggle_fullscreen(self):
        if self.fullscreen_mode:
            self.left_panel_ref.setVisible(True)
            self.btn_fullscreen.setText("🔍 Maximize 3D View")
            self.fullscreen_mode = False
        else:
            self.left_panel_ref.setVisible(False)
            self.btn_fullscreen.setText("◀ Show Controls")
            self.fullscreen_mode = True
        self.canvas.draw()
    
    def get_causal_pattern(self):
        patterns = {0: [True]*4, 1: [True, True, True, False], 
                   2: [True, True, True, False], 3: [True, True, False, False],
                   4: [True, False, True, False]}
        return patterns[self.combo_orientation.currentIndex()]
    
    def compute_configuration(self):
=======
        
        self.setLayout(main_layout)
        
        # Initial plot
        self.plot_empty_cone()
    
    def toggle_fullscreen(self):
        """Toggle between normal and fullscreen 3D view"""
        if self.fullscreen_mode:
            # Return to normal mode
            self.left_panel_ref.setVisible(True)
            self.splitter_ref.setStretchFactor(0, 1)
            self.splitter_ref.setStretchFactor(1, 3)
            self.btn_fullscreen.setText("🔍 Maximize 3D View")
            self.fullscreen_mode = False
        else:
            # Go to fullscreen mode
            self.left_panel_ref.setVisible(False)
            self.btn_fullscreen.setText("◀ Show Controls")
            self.fullscreen_mode = True
        
        self.canvas.draw()
    
    def get_causal_pattern(self):
        """Convert combo selection to causal pattern"""
        patterns = {
            0: [True, True, True, True],
            1: [True, True, True, False],
            2: [True, True, True, False],
            3: [True, True, False, False],
            4: [True, False, True, False]
        }
        return patterns[self.combo_orientation.currentIndex()]
    
    def compute_configuration(self):
        """Compute normals for current configuration"""
        # Stop previous thread if running
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.thread.wait()
        
        self.btn_compute.setEnabled(False)
        self.progress.setValue(0)
        
        params = {
            'j0': self.spin_j0.value(),
            'trials': self.spin_trials.value(),
            'causal_pattern': self.get_causal_pattern()
        }
        
        self.thread = ComputationThread("visualization", params)
<<<<<<< HEAD
        self.thread.progress.connect(lambda c, t: self.progress.setValue(c) or self.progress.setMaximum(t))
        self.thread.finished.connect(self.display_visualization)
        self.thread.error.connect(lambda e: self.info_text.setText(f"ERROR:\n{e}") or self.btn_compute.setEnabled(True))
        self.thread.start()
    
    def display_visualization(self, results):
=======
        self.thread.progress.connect(self.update_progress)
        self.thread.finished.connect(self.display_visualization)
        self.thread.error.connect(self.handle_error)
        self.thread.start()
    
    def update_progress(self, current, total):
        self.progress.setMaximum(total)
        self.progress.setValue(current)
    
    def handle_error(self, error_msg):
        self.info_text.setText(f"ERROR:\n{error_msg}")
        self.btn_compute.setEnabled(True)
        QMessageBox.warning(self, "Computation Error", 
                          "An error occurred during computation.")
    
    def display_visualization(self, results):
        """Display computed configuration"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        self.current_results = results
        self.update_plot()
        self.update_info()
        self.btn_compute.setEnabled(True)
    
    def update_info(self):
<<<<<<< HEAD
        if not self.current_results:
            return
        try:
            r = self.current_results
            pattern_str = "".join(["T+" if f else "T-" for f in r['causal_pattern']])
            text = f"Configuration: {pattern_str}\n"
            text += f"F = {r['F']:.3e}\n|det(G)| = {abs(r['det_G']):.3e}\n"
            closure = closure_constraint(r['normals'], r['spins'])
            text += f"||closure|| = {norm(closure):.3e}\n\nNormals:\n"
            for i, n in enumerate(r['normals']):
                text += f"n{i}: ({n[0]:+.3f}, {n[1]:+.3f}, {n[2]:+.3f}, {n[3]:+.3f})\n"
            self.info_text.setText(text)
        except Exception as e:
            self.info_text.setText(f"Error: {e}")
    
    def plot_empty_cone(self):
        try:
            ax = self.canvas.axes
            ax.clear()
            self.draw_light_cone(ax)
            ax.set_xlabel('x¹', fontsize=12, fontweight='bold', color='white')
            ax.set_ylabel('x²', fontsize=12, fontweight='bold', color='white')
            ax.set_zlabel('x⁰ (time)', fontsize=12, fontweight='bold', color='white')
            ax.set_title('Minkowski Light Cone', fontsize=12, fontweight='bold', color='white')
=======
        """Update info text"""
        if not self.current_results:
            return
        
        try:
            r = self.current_results
            pattern_str = "".join(["T+" if f else "T-" for f in r['causal_pattern']])
            
            text = f"Configuration: {pattern_str}\n"
            text += f"Objective F = {r['F']:.3e}\n"
            text += f"|det(G)| = {abs(r['det_G']):.3e}\n"
            
            # Check closure
            closure = closure_constraint(r['normals'], r['spins'])
            text += f"Closure norm = {norm(closure):.3e}\n"
            
            # Individual norms
            text += "\nNormal vectors (n⁰, n¹, n², n³):\n"
            for i, n in enumerate(r['normals']):
                text += f"n{i}: ({n[0]:+.3f}, {n[1]:+.3f}, {n[2]:+.3f}, {n[3]:+.3f})\n"
            
            self.info_text.setText(text)
        except Exception as e:
            self.info_text.setText(f"Error updating info: {str(e)}")
    
    def plot_empty_cone(self):
        """Plot empty light cone as initial state"""
        try:
            ax = self.canvas.axes
            ax.clear()
            
            self.draw_light_cone(ax)
            
            ax.set_xlabel('x¹', fontsize=12, fontweight='bold', color='white')
            ax.set_ylabel('x²', fontsize=12, fontweight='bold', color='white')
            ax.set_zlabel('x⁰ (time)', fontsize=12, fontweight='bold', color='white')
            ax.set_title('Minkowski Light Cone (Click "Compute Configuration" to start)', 
                        fontsize=12, fontweight='bold', color='white')
            
            # Style
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.grid(True, alpha=0.2)
            ax.tick_params(colors='white')
<<<<<<< HEAD
            self.canvas.draw()
        except Exception as e:
            logger.error(f"Error plotting cone: {e}")
    
    def draw_light_cone(self, ax, r_max=2.5, n_points=30):
=======
            
            self.canvas.draw()
        except Exception as e:
            print(f"Error plotting empty cone: {e}")
    
    def draw_light_cone(self, ax, r_max=2.5, n_points=30):
        """Draw the light cone surface"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        try:
            theta = np.linspace(0, 2*np.pi, n_points)
            r = np.linspace(0, r_max, n_points)
            Theta, R = np.meshgrid(theta, r)
<<<<<<< HEAD
=======
            
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            X = R * np.cos(Theta)
            Y = R * np.sin(Theta)
            Z_future = R
            Z_past = -R
<<<<<<< HEAD
            ax.plot_surface(X, Y, Z_future, alpha=0.15, color='cyan', linewidth=0)
            ax.plot_surface(X, Y, Z_past, alpha=0.15, color='magenta', linewidth=0)
        except Exception as e:
            logger.error(f"Error drawing cone: {e}")
    
    def draw_hyperboloid(self, ax, t_val=1.5, n_points=30):
=======
            
            # Future cone
            ax.plot_surface(X, Y, Z_future, alpha=0.15, color='cyan', 
                           linewidth=0, antialiased=True)
            
            # Past cone
            ax.plot_surface(X, Y, Z_past, alpha=0.15, color='magenta', 
                           linewidth=0, antialiased=True)
            
            # Cone edges
            for angle in np.linspace(0, 2*np.pi, 8, endpoint=False):
                x_line = r_max * np.cos(angle)
                y_line = r_max * np.sin(angle)
                ax.plot([0, x_line], [0, y_line], [0, r_max], 
                       color='cyan', alpha=0.4, linewidth=1)
                ax.plot([0, x_line], [0, y_line], [0, -r_max], 
                       color='magenta', alpha=0.4, linewidth=1)
        except Exception as e:
            print(f"Error drawing light cone: {e}")
    
    def draw_hyperboloid(self, ax, t_val=1.5, n_points=30):
        """Draw the timelike hyperboloid n·n = -1"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        try:
            theta = np.linspace(0, 2*np.pi, n_points)
            phi = np.linspace(0, np.pi, n_points//2)
            Theta, Phi = np.meshgrid(theta, phi)
<<<<<<< HEAD
=======
            
            # Future sheet
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            T_future = t_val
            R = np.sqrt(T_future**2 - 1)
            X = R * np.sin(Phi) * np.cos(Theta)
            Y = R * np.sin(Phi) * np.sin(Theta)
            Z = np.ones_like(X) * T_future
<<<<<<< HEAD
            ax.plot_surface(X, Y, Z, alpha=0.2, color='yellow', linewidth=0)
            Z_past = np.ones_like(X) * (-T_future)
            ax.plot_surface(X, Y, Z_past, alpha=0.2, color='orange', linewidth=0)
        except Exception as e:
            logger.error(f"Error drawing hyperboloid: {e}")
    
    def update_plot(self):
=======
            
            ax.plot_surface(X, Y, Z, alpha=0.2, color='yellow', 
                           linewidth=0, antialiased=True)
            
            # Past sheet
            T_past = -t_val
            Z_past = np.ones_like(X) * T_past
            ax.plot_surface(X, Y, Z_past, alpha=0.2, color='orange', 
                           linewidth=0, antialiased=True)
        except Exception as e:
            print(f"Error drawing hyperboloid: {e}")
    
    def update_plot(self):
        """Update 3D plot with current settings"""
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        try:
            ax = self.canvas.axes
            ax.clear()
            
<<<<<<< HEAD
            if self.check_cone.isChecked():
                self.draw_light_cone(ax)
            if self.check_hyperboloid.isChecked():
                self.draw_hyperboloid(ax)
            
            if self.current_results and self.check_normals.isChecked():
                for i, (n, future) in enumerate(zip(self.current_results['normals'], 
                                                    self.current_results['causal_pattern'])):
                    color = 'lime' if future else 'red'
                    scale = 1.5
                    ax.quiver(0, 0, 0, scale*n[1], scale*n[2], scale*n[0],
                             color=color, arrow_length_ratio=0.15, linewidth=2.5, alpha=0.9)
                    ax.text(scale*n[1], scale*n[2], scale*n[0], f'  n{i}', 
                           color=color, fontsize=10, fontweight='bold')
            
            if self.current_results and self.check_closure.isChecked():
                closure = closure_constraint(self.current_results['normals'], 
                                            self.current_results['spins'])
                scale = 0.1
                ax.quiver(0, 0, 0, scale*closure[1], scale*closure[2], scale*closure[0],
                         color='white', arrow_length_ratio=0.2, linewidth=3, alpha=0.8)
            
            ax.set_xlabel('x¹', fontsize=12, fontweight='bold', color='white')
            ax.set_ylabel('x²', fontsize=12, fontweight='bold', color='white')
            ax.set_zlabel('x⁰', fontsize=12, fontweight='bold', color='white')
            ax.set_xlim([-2.5, 2.5])
            ax.set_ylim([-2.5, 2.5])
            ax.set_zlim([-2.5, 2.5])
=======
            # Draw light cone
            if self.check_cone.isChecked():
                self.draw_light_cone(ax)
            
            # Draw hyperboloid
            if self.check_hyperboloid.isChecked():
                self.draw_hyperboloid(ax)
            
            # Draw normals if computed
            if self.current_results and self.check_normals.isChecked():
                normals = self.current_results['normals']
                spins = self.current_results['spins']
                pattern = self.current_results['causal_pattern']
                
                for i, (n, j, future) in enumerate(zip(normals, spins, pattern)):
                    # Extract spatial components for plotting
                    # We plot (n¹, n², n⁰) for better visualization
                    color = 'lime' if future else 'red'
                    label = f'n{i} (T+)' if future else f'n{i} (T-)'
                    
                    # Scale for visibility
                    scale = 1.5
                    ax.quiver(0, 0, 0, 
                             scale*n[1], scale*n[2], scale*n[0],
                             color=color, arrow_length_ratio=0.15, 
                             linewidth=2.5, alpha=0.9, label=label)
                    
                    # Add label
                    ax.text(scale*n[1], scale*n[2], scale*n[0], 
                           f'  n{i}', color=color, fontsize=10, 
                           fontweight='bold')
            
            # Draw closure vector
            if self.current_results and self.check_closure.isChecked():
                closure = closure_constraint(self.current_results['normals'], 
                                            self.current_results['spins'])
                
                # Scale for visibility
                scale = 0.1
                ax.quiver(0, 0, 0,
                         scale*closure[1], scale*closure[2], scale*closure[0],
                         color='white', arrow_length_ratio=0.2,
                         linewidth=3, alpha=0.8, label='Σ jᵢnᵢ (closure)')
                
                ax.text(scale*closure[1], scale*closure[2], scale*closure[0],
                       '  Closure', color='white', fontsize=10, fontweight='bold')
            
            # Labels and title
            ax.set_xlabel('x¹', fontsize=12, fontweight='bold', color='white')
            ax.set_ylabel('x²', fontsize=12, fontweight='bold', color='white')
            ax.set_zlabel('x⁰ (time)', fontsize=12, fontweight='bold', color='white')
            
            if self.current_results:
                pattern_str = "".join(["T+" if f else "T-" 
                                      for f in self.current_results['causal_pattern']])
                title = f'Configuration: {pattern_str}'
            else:
                title = 'Minkowski Light Cone'
            
            ax.set_title(title, fontsize=12, fontweight='bold', color='white')
            
            # Set limits
            ax.set_xlim([-2.5, 2.5])
            ax.set_ylim([-2.5, 2.5])
            ax.set_zlim([-2.5, 2.5])
            
            # Style
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.grid(True, alpha=0.2)
            ax.tick_params(colors='white')
            
<<<<<<< HEAD
            self.update_view()
            self.canvas.draw()
        except Exception as e:
            logger.error(f"Error updating plot: {e}")
    
    def update_view(self):
        try:
            elev = self.slider_elev.value()
            azim = self.slider_azim.value()
            self.label_elev.setText(f"{elev}°")
            self.label_azim.setText(f"{azim}°")
            self.canvas.axes.view_init(elev=elev, azim=azim)
            self.canvas.draw()
        except:
            pass

# =============================================================================
# THEORY TAB
# =============================================================================

class TheoryTab(QWidget):
    """Tab showing theoretical background"""
    
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
=======
            # Legend
            if self.current_results and self.check_normals.isChecked():
                ax.legend(loc='upper left', fontsize=8, facecolor='#2b2b2b', 
                         edgecolor='white', labelcolor='white')
            
            # Update view
            self.update_view()
            
            self.canvas.draw()
        except Exception as e:
            print(f"Error updating plot: {e}")
    
    def update_view(self):
        """Update 3D view angle"""
        try:
            elev = self.slider_elev.value()
            azim = self.slider_azim.value()
            
            # Update labels
            self.label_elev.setText(f"{elev}°")
            self.label_azim.setText(f"{azim}°")
            
            self.canvas.axes.view_init(elev=elev, azim=azim)
            self.canvas.draw()
        except Exception as e:
            print(f"Error updating view: {e}")

class TheoryTab(QWidget):
    """Tab showing theoretical background"""
    def __init__(self):
        super().__init__()
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        text = QTextEdit()
        text.setReadOnly(True)
        
        theory_text = """
╔════════════════════════════════════════════════════════════════════════════╗
║                    THEORETICAL BACKGROUND                                  ║
║  Causally Confined Euclidean Saddles in Spin-Foam Quantum Gravity         ║
╚════════════════════════════════════════════════════════════════════════════╝

<<<<<<< HEAD
VERSION 5.0 - PRODUCTION READY

MAJOR IMPROVEMENTS:
✅ Centralized configuration (no magic numbers)
✅ Robust input validation
✅ Timeout protection for optimization
✅ Memory-efficient sensitivity analysis
✅ Export results (JSON/CSV)
✅ Complete documentation
✅ Structured logging
✅ Reduced code duplication

1. MAIN RESULT (Section 2.3.3)
──────────────────────────────────────────────────────────────────────────────
    S_E = ℏj C({α})  →  A_v ~ exp(-j C({α}))

Euclidean action scales linearly with spin, leading to exponential suppression.

2. GEOMETRIC OBSTRUCTION (Appendix A)
──────────────────────────────────────────────────────────────────────────────
Proposition A.4: For MIXED temporal orientations, NO real Lorentzian solution.

Numerical Verification:
- Uniform: det(G) ~ 10^{-60}  ✓ (closure achieved)
- Mixed:   det(G) ~ 10^{-15}  ✗ (closure fails)

3. 3D VISUALIZATION
──────────────────────────────────────────────────────────────────────────────
- Future cone (cyan): x⁰ > 0
- Past cone (magenta): x⁰ < 0
- Hyperboloid: n·n = -1
- Normal vectors with color-coded orientations

The visualization makes the causal obstruction VISIBLE.

4. EXPORT FUNCTIONALITY (NEW IN v5.0)
──────────────────────────────────────────────────────────────────────────────
All analysis tabs now support exporting results to JSON and CSV formats.

5. REFERENCES
──────────────────────────────────────────────────────────────────────────────
[1] C. Rovelli, Quantum Gravity (Cambridge, 2004)
[2] Barrett et al., Class. Quantum Grav. 26, 165009 (2009)
[3] Han, Krajewski, Rovelli, CQG 30, 165012 (2013)

See paper for complete bibliography.

╔════════════════════════════════════════════════════════════════════════════╗
║  "The saddle emerges as a necessity, not a postulate." — Section A.7      ║
╚════════════════════════════════════════════════════════════════════════════╝
=======
1. CONTEXT
──────────────────────────────────────────────────────────────────────────────
This application provides numerical verification for Appendix A of the paper:

    "Causally Confined Euclidean Saddles in Spin-Foam Quantum Gravity"
    Mário Sérgio Guilherme Junior (2025)

The paper investigates the microscopic origin of exponential suppression weights
associated with localized quantum-geometric processes in covariant loop quantum
gravity.

2. MAIN RESULT
──────────────────────────────────────────────────────────────────────────────
Key equation from Section 2.3.3:

    S_E = ℏj C({α})     where C({α}) = Σ α_ab χ_ab({α})
                                      a<b

This shows that the Euclidean action scales linearly with spin j, leading to
exponential suppression:

    A_v ~ exp(-S_E/ℏ) = exp(-j C({α}))

3. GEOMETRIC OBSTRUCTION (Appendix A)
──────────────────────────────────────────────────────────────────────────────
Proposition A.4: For boundary data with MIXED temporal orientations, and for
generic spins (outside a measure-zero subset), there exists NO real Lorentzian
solution to the closure and gluing equations.

Analytical Argument:
• Future-pointing timelike vectors form a convex cone in ℝ^{3,1}
• Sum of future-pointing vectors is necessarily future-pointing
• Cannot cancel with past-pointing vector except in degenerate cases

Numerical Verification:
• Closure: Σ j_i n_i = 0  requires linear dependence → det(G) ≈ 0
• For uniform orientation: det(G) ~ 10^{-60}  ✓ (closure achieved)
• For mixed orientation:   det(G) ~ 10^{-15}  ✗ (closure fails)

4. 3D VISUALIZATION (v4.1)
──────────────────────────────────────────────────────────────────────────────
The Minkowski light cone visualization shows:

• Future cone (cyan): {(x⁰,x¹,x²,x³) : (x⁰)² = (x¹)² + (x²)² + (x³)², x⁰>0}
• Past cone (magenta): {(x⁰,x¹,x²,x³) : (x⁰)² = (x¹)² + (x²)² + (x³)², x⁰<0}
• Hyperboloid: {n : n·n = -1} (timelike unit vectors)
• Normal vectors: n_i satisfying closure Σ j_i n_i = 0

Color coding:
• Green (lime) arrows: Future-pointing normals (T+)
• Red arrows: Past-pointing normals (T-)
• White arrow: Closure vector (should be ~0 for valid geometries)

The visualization makes the causal obstruction VISIBLE:
→ Uniform (T+): All green arrows can close (lie in same cone)
→ Mixed (T+/T-): Green + red arrows cannot close (opposite cones)

5. PHYSICAL INTERPRETATION
──────────────────────────────────────────────────────────────────────────────
When no real Lorentzian geometry exists:
• Path integral requires analytic continuation: θ → iχ
• Action becomes purely imaginary: S_Regge → iS_E
• Amplitude becomes non-oscillatory: exp(iS/ℏ) → exp(-S_E/ℏ)
• This describes a TUNNELING process in the space of discrete geometries

Connection to MQT Framework:
The exponential weight ε exp(-S_E/ℏ) in the Momentary Quantum Tunneling
framework emerges directly from this spin-foam mechanism.

6. SCOPE AND LIMITATIONS
──────────────────────────────────────────────────────────────────────────────
This is a PROOF-OF-CONCEPT demonstration:

✓ Shows mechanism is viable in toy model
✓ Provides numerical evidence for analytical obstruction
✓ Demonstrates semiclassical control (S_E ~ j)
✓ Visualizes geometric structure of obstruction

Open questions:
• Universality beyond toy model
• Triangulation independence
• Full spin-foam sum
• Microscopic derivation of prefactor ε

7. REFERENCES
──────────────────────────────────────────────────────────────────────────────
[1] C. Rovelli, Quantum Gravity (Cambridge, 2004)
[2] Barrett et al., "EPRL amplitude" (Class. Quantum Grav. 26, 165009, 2009)
[3] Han, Krajewski, Rovelli, "Complex saddle points" (CQG 30, 165012, 2013)
[4] Coleman, "Fate of false vacuum" (Phys. Rev. D 15, 2929, 1977)

See paper for complete bibliography and detailed derivations.

╔════════════════════════════════════════════════════════════════════════════╗
║  "The saddle emerges as a necessity, not a postulate."                     ║
║                                                    — Section A.7            ║
╚════════════════════════════════════════════════════════════════════════════╝

VERSION 4.2 UPDATES:
──────────────────────────────────────────────────────────────────────────────
✓ Fixed thread safety issues
✓ Added proper error handling
✓ Improved progress tracking
✓ Added stop button for long computations
✓ Better memory management
✓ More robust optimization
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        """
        
        text.setPlainText(theory_text)
        text.setFont(QFont("Courier", 9))
<<<<<<< HEAD
=======
        
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        layout.addWidget(text)
        self.setLayout(layout)

# =============================================================================
# MAIN WINDOW
# =============================================================================

class MainWindow(QMainWindow):
    """Main application window"""
<<<<<<< HEAD
    
=======
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
    def __init__(self):
        super().__init__()
        self.init_ui()
        
    def init_ui(self):
<<<<<<< HEAD
        self.setWindowTitle("Spin-Foam Analysis v5.0 PRODUCTION")
        self.setGeometry(50, 50, Config.MAIN_WINDOW_WIDTH, Config.MAIN_WINDOW_HEIGHT)
        self.apply_style()
        
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout()
        
=======
        self.setWindowTitle("Spin-Foam Euclidean Saddles Analysis v4.3 [PROGRESS FIX]")
        self.setGeometry(50, 50, 1600, 1000)  # Larger window size
        
        # Apply dark theme
        self.apply_style()
        
        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        
        layout = QVBoxLayout()
        
        # Header
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        header = QLabel("Causally Confined Euclidean Saddles in Spin-Foam Quantum Gravity")
        header.setFont(QFont("Arial", 16, QFont.Bold))
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
<<<<<<< HEAD
        subheader = QLabel("v5.0 PRODUCTION • M. S. Guilherme Junior")
=======
        subheader = QLabel("Numerical Verification Tool • Version 4.3 • M. S. Guilherme Junior")
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        subheader.setFont(QFont("Arial", 10))
        subheader.setAlignment(Qt.AlignCenter)
        layout.addWidget(subheader)
        
<<<<<<< HEAD
        version_info = QLabel("✅ Config ✅ Validation ✅ Timeout ✅ Export ✅ Logging ✅ Docs")
=======
        version_info = QLabel("✅ PROGRESS FIX: Bars now reach 100% correctly")
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        version_info.setFont(QFont("Arial", 9))
        version_info.setAlignment(Qt.AlignCenter)
        version_info.setStyleSheet("color: #00ff00;")
        layout.addWidget(version_info)
        
<<<<<<< HEAD
        tabs = QTabWidget()
        tabs.addTab(TheoryTab(), "📚 Theory")
        tabs.addTab(BasicCheckTab(), "✓ Basic Check")
        tabs.addTab(ScalingTab(), "📈 Scaling")
        tabs.addTab(SensitivityTab(), "🎯 Sensitivity")
        tabs.addTab(VisualizationTab(), "🔮 3D Viz")
        layout.addWidget(tabs)
        
        footer = QLabel("Production-Ready • Full Documentation • Robust Error Handling")
=======
        # Tabs
        tabs = QTabWidget()
        tabs.addTab(TheoryTab(), "📚 Theory")
        tabs.addTab(BasicCheckTab(), "✓ Basic Check (Table 1)")
        tabs.addTab(ScalingTab(), "📈 Scaling Analysis")
        tabs.addTab(SensitivityTab(), "🎯 Sensitivity")
        tabs.addTab(VisualizationTab(), "🔮 3D Visualization")
        
        layout.addWidget(tabs)
        
        # Footer
        footer = QLabel("Proof-of-Concept • See paper for complete theoretical framework")
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        footer.setFont(QFont("Arial", 9))
        footer.setAlignment(Qt.AlignCenter)
        layout.addWidget(footer)
        
        central.setLayout(layout)
    
    def apply_style(self):
<<<<<<< HEAD
        """Apply dark theme stylesheet"""
        style = """
        QMainWindow { background-color: #2b2b2b; }
        QLabel { color: #ffffff; }
=======
        """Apply custom stylesheet"""
        style = """
        QMainWindow {
            background-color: #2b2b2b;
        }
        QLabel {
            color: #ffffff;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        QGroupBox {
            color: #ffffff;
            border: 2px solid #555555;
            border-radius: 5px;
            margin-top: 10px;
            font-weight: bold;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
<<<<<<< HEAD
            padding: 0 5px;
=======
            padding: 0 5px 0 5px;
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        }
        QPushButton {
            background-color: #0d7377;
            color: white;
            border: none;
            padding: 8px;
            border-radius: 4px;
            font-weight: bold;
        }
<<<<<<< HEAD
        QPushButton:hover { background-color: #14a085; }
        QPushButton:pressed { background-color: #0a5f63; }
=======
        QPushButton:hover {
            background-color: #14a085;
        }
        QPushButton:pressed {
            background-color: #0a5f63;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        QPushButton:disabled {
            background-color: #555555;
            color: #888888;
        }
        QTextEdit {
            background-color: #1e1e1e;
            color: #00ff00;
            border: 1px solid #555555;
            border-radius: 3px;
        }
        QSpinBox, QDoubleSpinBox {
            background-color: #3c3c3c;
            color: #ffffff;
            border: 1px solid #555555;
            border-radius: 3px;
            padding: 3px;
        }
        QComboBox {
            background-color: #3c3c3c;
            color: #ffffff;
            border: 1px solid #555555;
            border-radius: 3px;
            padding: 3px;
        }
<<<<<<< HEAD
=======
        QComboBox::drop-down {
            border: none;
        }
        QComboBox::down-arrow {
            image: none;
            border-left: 5px solid transparent;
            border-right: 5px solid transparent;
            border-top: 5px solid #ffffff;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        QProgressBar {
            border: 1px solid #555555;
            border-radius: 3px;
            text-align: center;
            color: white;
        }
<<<<<<< HEAD
        QProgressBar::chunk { background-color: #0d7377; }
=======
        QProgressBar::chunk {
            background-color: #0d7377;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        QTabWidget::pane {
            border: 1px solid #555555;
            background-color: #2b2b2b;
        }
        QTabBar::tab {
            background-color: #3c3c3c;
            color: #ffffff;
            padding: 8px 16px;
            margin-right: 2px;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }
<<<<<<< HEAD
        QTabBar::tab:selected { background-color: #0d7377; }
        QTabBar::tab:hover { background-color: #4a4a4a; }
        QCheckBox { color: #ffffff; }
=======
        QTabBar::tab:selected {
            background-color: #0d7377;
        }
        QTabBar::tab:hover {
            background-color: #4a4a4a;
        }
        QCheckBox {
            color: #ffffff;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        QCheckBox::indicator {
            width: 18px;
            height: 18px;
            border: 2px solid #555555;
            border-radius: 3px;
            background-color: #3c3c3c;
        }
<<<<<<< HEAD
        QCheckBox::indicator:checked { background-color: #0d7377; }
=======
        QCheckBox::indicator:checked {
            background-color: #0d7377;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        QSlider::groove:horizontal {
            height: 8px;
            background: #3c3c3c;
            border-radius: 4px;
        }
        QSlider::handle:horizontal {
            background: #0d7377;
            width: 18px;
            margin: -5px 0;
            border-radius: 9px;
        }
<<<<<<< HEAD
        QSlider::handle:horizontal:hover { background: #14a085; }
=======
        QSlider::handle:horizontal:hover {
            background: #14a085;
        }
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        """
        self.setStyleSheet(style)
    
    def closeEvent(self, event):
<<<<<<< HEAD
        """Ensure all threads stop on close"""
        logger.info("Application closing...")
=======
        """Handle window close event"""
        # Ensure all threads are stopped
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
        for widget in self.findChildren(QWidget):
            if hasattr(widget, 'thread') and widget.thread:
                if widget.thread.isRunning():
                    widget.thread.stop()
<<<<<<< HEAD
                    widget.thread.wait(Config.THREAD_WAIT_TIMEOUT_MS)
        event.accept()
        logger.info("Application closed successfully")
=======
                    widget.thread.wait(1000)  # Wait up to 1 second
        event.accept()
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d

# =============================================================================
# MAIN
# =============================================================================

def main():
<<<<<<< HEAD
    """Main application entry point"""
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    logger.info("Launching main window...")
    window = MainWindow()
    window.show()
    
    logger.info("Application ready")
    exit_code = app.exec_()
    
    logger.info(f"Application exited with code {exit_code}")
    logger.info("="*80)
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
=======
    app = QApplication(sys.argv)
    app.setStyle('Fusion')  # Modern look
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
>>>>>>> 4eb069c8c0b16556ebc43259d49bd7f37ef6443d
