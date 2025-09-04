import ROOT
import numpy as np
import sys
import os
import argparse
from array import array

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes

from common_functions import *
from binning_and_selections import *

def Proj2D_Y(h,xmin,xmax,hname="XXX"):

    # project 2D histogram into 1D along Y

    imin=h.GetXaxis().FindBin(xmin)
    imax=h.GetXaxis().FindBin(xmax)-1
    
    proj_y=h.ProjectionY(hname, imin, imax)
    ROOT.SetOwnership(proj_y,True)

    return proj_y

def Proj2D_X(h,ymin,ymax,hname="XXX",Debug=False):

    # project 2D histogram into 1D along Y

    imin=h.GetYaxis().FindBin(ymin)
    imax=h.GetYaxis().FindBin(ymax)-1

    proj_x=h.ProjectionX(hname, imin, imax)
    ROOT.SetOwnership(proj_x,True)

    return proj_x

def flatten_2d_histogram_with_overflow(hist_2d, name_suffix="_flattened", include_overflow=True):
    """
    Convert 2D histogram to 1D by flattening, with option to include overflow/underflow.
    """
    
    if include_overflow:
        # Include overflow/underflow bins
        nx = hist_2d.GetNbinsX() + 2  # +2 for under/overflow
        ny = hist_2d.GetNbinsY() + 2  # +2 for under/overflow
        i_start, i_end = 0, nx
        j_start, j_end = 0, ny
        print(f"Including overflow: {nx}x{ny} bins (with under/overflow)")
    else:
        # Normal bins only
        nx = hist_2d.GetNbinsX()
        ny = hist_2d.GetNbinsY()
        i_start, i_end = 1, nx + 1
        j_start, j_end = 1, ny + 1
        print(f"Normal bins only: {nx}x{ny} bins")
    
    total_bins = nx * ny
    
    print(f"Flattening 2D histogram -> 1D with {total_bins} bins")
    
    # Create 1D histogram
    flat_name = hist_2d.GetName() + name_suffix
    hist_1d = ROOT.TH1D(flat_name, hist_2d.GetTitle(), total_bins, 0.5, total_bins + 0.5)
    
    # Extract bin edges (for normal bins only, overflow doesn't have meaningful edges)
    x_edges = []
    y_edges = []
    
    if include_overflow:
        # For overflow bins, we need to handle edges carefully
        # We need to create extended edges that include the original nominal edges
        # plus extensions for underflow/overflow
        
        # Get nominal bin edges first
        x_nominal_edges = []
        y_nominal_edges = []
        
        for i in range(hist_2d.GetNbinsX() + 1):
            x_nominal_edges.append(hist_2d.GetXaxis().GetBinLowEdge(i + 1))
        
        for j in range(hist_2d.GetNbinsY() + 1):
            y_nominal_edges.append(hist_2d.GetYaxis().GetBinLowEdge(j + 1))
        
        # Create extended edges: [underflow_edge] + [nominal_edges] + [overflow_edge]
        # Use a reasonable extension based on bin width
        x_first_width = x_nominal_edges[1] - x_nominal_edges[0]
        x_last_width = x_nominal_edges[-1] - x_nominal_edges[-2]
        y_first_width = y_nominal_edges[1] - y_nominal_edges[0]
        y_last_width = y_nominal_edges[-1] - y_nominal_edges[-2]
        
        # Build extended edge arrays
        x_edges = [x_nominal_edges[0] - x_first_width] + x_nominal_edges + [x_nominal_edges[-1] + x_last_width]
        y_edges = [y_nominal_edges[0] - y_first_width] + y_nominal_edges + [y_nominal_edges[-1] + y_last_width]
        
        print(f"Extended X edges: {len(x_edges)} total ({len(x_nominal_edges)} nominal + 2 overflow)")
        print(f"Extended Y edges: {len(y_edges)} total ({len(y_nominal_edges)} nominal + 2 overflow)")
    else:
        # Normal bin edges
        for i in range(nx + 1):
            x_edges.append(hist_2d.GetXaxis().GetBinLowEdge(i + 1))
        for j in range(ny + 1):
            y_edges.append(hist_2d.GetYaxis().GetBinLowEdge(j + 1))
    
    # Enhanced mapping information
    mapping = {
        'nx': nx,
        'ny': ny,
        'include_overflow': include_overflow,
        'x_edges': x_edges,
        'y_edges': y_edges,
        'x_min': hist_2d.GetXaxis().GetXmin(),
        'x_max': hist_2d.GetXaxis().GetXmax(),
        'y_min': hist_2d.GetYaxis().GetXmin(),
        'y_max': hist_2d.GetYaxis().GetXmax()
    }
    
    # Fill 1D histogram
    for i in range(i_start, i_end):
        for j in range(j_start, j_end):
            content = hist_2d.GetBinContent(i, j)
            error = hist_2d.GetBinError(i, j)
            
            # Calculate flat bin index
            if include_overflow:
                flat_bin = i * ny + j + 1  # i and j start from 0
            else:
                flat_bin = (i - 1) * ny + j  # i and j start from 1
            
            hist_1d.SetBinContent(flat_bin, content)
            hist_1d.SetBinError(flat_bin, error)
    
    print(f"Original integral: {hist_2d.Integral():.6f}")
    print(f"Flattened integral: {hist_1d.Integral():.6f}")
    
    if include_overflow:
        # Check overflow/underflow content
        overflow_content = 0.0
        for i in [0, nx-1]:  # First and last i (underflow/overflow)
            for j in range(ny):
                overflow_content += hist_1d.GetBinContent(i * ny + j + 1)
        for j in [0, ny-1]:  # First and last j (underflow/overflow) 
            for i in range(1, nx-1):  # Avoid double counting corners
                overflow_content += hist_1d.GetBinContent(i * ny + j + 1)
        
        print(f"Overflow/underflow content: {overflow_content:.6f}")
    
    return hist_1d, mapping

def flatten_2d_histogram(hist_2d, name_suffix="_flattened"):
    """
    Convert 2D histogram to 1D by flattening (row-major order).
    """
    
    nx = hist_2d.GetNbinsX()
    ny = hist_2d.GetNbinsY()
    total_bins = nx * ny
    
    print(f"Flattening 2D histogram {nx}x{ny} -> 1D with {total_bins} bins")
    
    # Create 1D histogram
    flat_name = hist_2d.GetName() + name_suffix
    hist_1d = ROOT.TH1D(flat_name, hist_2d.GetTitle(), total_bins, 0.5, total_bins + 0.5)
    
    # FIXED: Extract actual bin edges for perfect reconstruction
    x_edges = []
    y_edges = []
    
    # Get all X bin edges (including upper edge of last bin)
    for i in range(nx + 1):
        x_edges.append(hist_2d.GetXaxis().GetBinLowEdge(i + 1))
    
    # Get all Y bin edges (including upper edge of last bin)
    for j in range(ny + 1):
        y_edges.append(hist_2d.GetYaxis().GetBinLowEdge(j + 1))
    
    # Enhanced mapping information with actual bin edges
    mapping = {
        'nx': nx,
        'ny': ny,
        'x_edges': x_edges,  # Actual bin edges for X-axis
        'y_edges': y_edges,  # Actual bin edges for Y-axis
        'x_min': hist_2d.GetXaxis().GetXmin(),
        'x_max': hist_2d.GetXaxis().GetXmax(),
        'y_min': hist_2d.GetYaxis().GetXmin(),
        'y_max': hist_2d.GetYaxis().GetXmax()
    }
    
    # Fill 1D histogram (row-major: i + j*nx + 1)
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            content = hist_2d.GetBinContent(i, j)
            error = hist_2d.GetBinError(i, j)
            
            #flat_bin = (j - 1) * nx + i
            flat_bin = (i - 1) * ny + j
            
            hist_1d.SetBinContent(flat_bin, content)
            hist_1d.SetBinError(flat_bin, error)
    
    print(f"Created flattened histogram: {hist_1d.GetNbinsX()} bins")
    print(f"Original integral: {hist_2d.Integral():.6f}")
    print(f"Flattened integral: {hist_1d.Integral():.6f}")
    print(f"Preserved X bin edges: {len(x_edges)} edges")
    print(f"Preserved Y bin edges: {len(y_edges)} edges")
    
    return hist_1d, mapping

def unflatten_1d_to_2d_nominal_only(hist_1d, mapping, name_suffix="_unflattened"):
    """
    Convert flattened 1D histogram back to 2D, keeping ONLY nominal range bins.
    Overflow/underflow bins are captured in the response matrix but ignored in final result.
    """
    
    nx = mapping['nx']
    ny = mapping['ny']
    include_overflow = mapping.get('include_overflow', False)
    
    if include_overflow:
        # Original histogram dimensions (without overflow)
        nx_nominal = nx - 2
        ny_nominal = ny - 2
        print(f"Unflattening with overflow -> nominal 2D {nx_nominal}x{ny_nominal} (ignoring overflow bins)")
    else:
        nx_nominal = nx
        ny_nominal = ny
        print(f"Unflattening 1D histogram -> 2D {nx_nominal}x{ny_nominal}")
    
    # Create 2D histogram with ONLY nominal range
    unflatten_name = hist_1d.GetName() + name_suffix
    
    if 'x_edges' in mapping and 'y_edges' in mapping:
        x_edges_array = array('d', mapping['x_edges'])
        y_edges_array = array('d', mapping['y_edges'])
        
        if include_overflow:
            # FIXED: Correctly extract only the nominal bin edges
            # The overflow array has structure: [underflow_edge, nominal_edges..., overflow_edge]
            # For nx_nominal bins, we need nx_nominal+1 edges
            # Skip first and last edges which are overflow extensions
            
            if len(x_edges_array) != nx + 1 or len(y_edges_array) != ny + 1:
                print(f"WARNING: Edge array size mismatch. X: {len(x_edges_array)} vs {nx+1}, Y: {len(y_edges_array)} vs {ny+1}")
            
            # Extract nominal edges: skip first (underflow) and last (overflow)
            x_nominal_edges = array('d', x_edges_array[1:nx_nominal+2])  # [1:nx-1] = nominal edges
            y_nominal_edges = array('d', y_edges_array[1:ny_nominal+2])  # [1:ny-1] = nominal edges
            
            print(f"Extracted {len(x_nominal_edges)-1} nominal X bins from {len(x_edges_array)-1} total")
            print(f"Extracted {len(y_nominal_edges)-1} nominal Y bins from {len(y_edges_array)-1} total")
            print(f"X range: [{x_nominal_edges[0]:.6f}, {x_nominal_edges[-1]:.6f}]")
            print(f"Y range: [{y_nominal_edges[0]:.6f}, {y_nominal_edges[-1]:.6f}]")
            
            hist_2d = ROOT.TH2D(
                unflatten_name, hist_1d.GetTitle(),
                nx_nominal, x_nominal_edges,
                ny_nominal, y_nominal_edges
            )
        else:
            hist_2d = ROOT.TH2D(
                unflatten_name, hist_1d.GetTitle(),
                nx_nominal, x_edges_array,
                ny_nominal, y_edges_array
            )
    else:
        # Fallback to uniform binning - always use nominal range
        hist_2d = ROOT.TH2D(
            unflatten_name, hist_1d.GetTitle(),
            nx_nominal, mapping['x_min'], mapping['x_max'],
            ny_nominal, mapping['y_min'], mapping['y_max']
        )
    
    # Fill 2D histogram - ONLY nominal bins
    nominal_content = 0.0
    overflow_content = 0.0
    
    for flat_bin in range(1, hist_1d.GetNbinsX() + 1):
        content = hist_1d.GetBinContent(flat_bin)
        error = hist_1d.GetBinError(flat_bin)
        
        # Reverse the flattening
        linear_index = flat_bin - 1
        
        if include_overflow:
            i = linear_index // ny  # 0-based, includes overflow
            j = linear_index % ny   # 0-based, includes overflow
            
            # Only fill nominal bins (exclude overflow: i=0, i=nx-1, j=0, j=ny-1)
            if i > 0 and i < nx-1 and j > 0 and j < ny-1:
                # Convert to 1-based nominal histogram coordinates
                # Subtract 1 because we're skipping the underflow bin
                i_nominal = i         # i=1..nx-2 maps to 1..nx_nominal  
                j_nominal = j         # j=1..ny-2 maps to 1..ny_nominal
                hist_2d.SetBinContent(i_nominal, j_nominal, content)
                hist_2d.SetBinError(i_nominal, j_nominal, error)
                nominal_content += content
            else:
                # This is overflow/underflow content - just track it
                overflow_content += content
        else:
            # No overflow - straightforward mapping
            i = (linear_index // ny) + 1  # 1-based
            j = (linear_index % ny) + 1   # 1-based
            hist_2d.SetBinContent(i, j, content)
            hist_2d.SetBinError(i, j, error)
            nominal_content += content
    
    print(f"Nominal range integral: {nominal_content:.6f}")
    if include_overflow and overflow_content > 0:
        print(f"Overflow/underflow content ignored: {overflow_content:.6f}")
        print(f"Fraction kept in nominal range: {nominal_content/(nominal_content + overflow_content):.3f}")
    
    print(f"Final 2D histogram integral: {hist_2d.Integral():.6f}")
    
    return hist_2d

def unflatten_1d_to_2d(hist_1d, mapping, name_suffix="_unflattened"):
    """
    Convert flattened 1D histogram back to 2D.
    """
    
    nx = mapping['nx']
    ny = mapping['ny']
    
    print(f"Unflattening 1D histogram -> 2D {nx}x{ny}")
    
    # Create 2D histogram with exact bin edges
    unflatten_name = hist_1d.GetName() + name_suffix
    
    # FIXED: Use actual bin edges if available
    if 'x_edges' in mapping and 'y_edges' in mapping:
        print("Using preserved bin edges for non-uniform binning")
        
        x_edges_array = array('d', mapping['x_edges'])
        y_edges_array = array('d', mapping['y_edges'])
        
        hist_2d = ROOT.TH2D(
            unflatten_name, hist_1d.GetTitle(),
            nx, x_edges_array,  # Use exact X bin edges
            ny, y_edges_array   # Use exact Y bin edges
        )
        
        # Verify bin edges were set correctly
        print(f"X-axis: {len(mapping['x_edges'])} edges, range [{x_edges_array[0]:.6f}, {x_edges_array[-1]:.6f}]")
        print(f"Y-axis: {len(mapping['y_edges'])} edges, range [{y_edges_array[0]:.6f}, {y_edges_array[-1]:.6f}]")
        
    else:
        print("WARNING: Using uniform binning (bin edges not preserved)")
        hist_2d = ROOT.TH2D(
            unflatten_name, hist_1d.GetTitle(),
            nx, mapping['x_min'], mapping['x_max'],
            ny, mapping['y_min'], mapping['y_max']
        )
    
    # Fill 2D histogram
    for flat_bin in range(1, hist_1d.GetNbinsX() + 1):
        content = hist_1d.GetBinContent(flat_bin)
        error = hist_1d.GetBinError(flat_bin)
        
        # Reverse the flattening: column-major
        linear_index = flat_bin - 1
        #i = (linear_index % nx) + 1
        #j = (linear_index // nx) + 1

        i = (linear_index // ny) + 1
        j = (linear_index % ny) + 1
        
        hist_2d.SetBinContent(i, j, content)
        hist_2d.SetBinError(i, j, error)
    
    print(f"Unflattened integral: {hist_2d.Integral():.6f}")
    
    return hist_2d

def create_2d_response_matrix_from_thnf_with_overflow(thnf_4d, name_suffix="_response_2d", include_overflow=True):
    """
    Convert 4D THnF to 2D TH2 response matrix, optionally including overflow/underflow.
    """
    
    print(f"Converting 4D THnF to 2D response matrix (include_overflow={include_overflow})...")
    
    if thnf_4d.GetNdimensions() != 4:
        raise ValueError(f"Expected 4D THnF, got {thnf_4d.GetNdimensions()}D")
    
    # Get dimensions
    axis_reco_x = thnf_4d.GetAxis(0)
    axis_reco_y = thnf_4d.GetAxis(1)
    axis_gen_x = thnf_4d.GetAxis(2)
    axis_gen_y = thnf_4d.GetAxis(3)
    
    if include_overflow:
        # Include overflow/underflow bins
        nx_reco = axis_reco_x.GetNbins() + 2
        ny_reco = axis_reco_y.GetNbins() + 2
        nx_gen = axis_gen_x.GetNbins() + 2
        ny_gen = axis_gen_y.GetNbins() + 2
        
        reco_start, reco_end = 0, nx_reco
        gen_start, gen_end = 0, nx_gen
        reco_j_start, reco_j_end = 0, ny_reco
        gen_j_start, gen_j_end = 0, ny_gen
    else:
        # Normal bins only
        nx_reco = axis_reco_x.GetNbins()
        ny_reco = axis_reco_y.GetNbins()
        nx_gen = axis_gen_x.GetNbins()
        ny_gen = axis_gen_y.GetNbins()
        
        reco_start, reco_end = 1, nx_reco + 1
        gen_start, gen_end = 1, nx_gen + 1
        reco_j_start, reco_j_end = 1, ny_reco + 1
        gen_j_start, gen_j_end = 1, ny_gen + 1
    
    total_reco_bins = nx_reco * ny_reco
    total_gen_bins = nx_gen * ny_gen
    
    print(f"THnF dimensions: reco({nx_reco}x{ny_reco}) -> gen({nx_gen}x{ny_gen})")
    print(f"Flattened response matrix: {total_reco_bins} x {total_gen_bins}")
    
    # Create 2D response matrix
    response_name = thnf_4d.GetName() + name_suffix
    response_2d = ROOT.TH2D(
        response_name, "2D Response Matrix",
        total_reco_bins, 0.5, total_reco_bins + 0.5,  # Reco (X-axis)
        total_gen_bins, 0.5, total_gen_bins + 0.5      # Gen (Y-axis)
    )
    
    print("Filling 2D response matrix...")
    
    filled_bins = 0
    coords = array('i', [0, 0, 0, 0])
    total_content = 0.0
    overflow_content = 0.0
    
    # Iterate through all bins (including overflow if requested)
    for i_reco in range(reco_start, reco_end):
        for j_reco in range(reco_j_start, reco_j_end):
            # Calculate flat reco bin
            if include_overflow:
                flat_reco_bin = i_reco * ny_reco + j_reco + 1
            else:
                flat_reco_bin = (i_reco - 1) * ny_reco + j_reco
            
            for i_gen in range(gen_start, gen_end):
                for j_gen in range(gen_j_start, gen_j_end):
                    # Calculate flat gen bin
                    if include_overflow:
                        flat_gen_bin = i_gen * ny_gen + j_gen + 1
                    else:
                        flat_gen_bin = (i_gen - 1) * ny_gen + j_gen
                    
                    # For THnF, we need to handle overflow bins differently
                    if include_overflow:
                        # Map 0-based indices to THnF overflow convention
                        if i_reco == 0:
                            coords[0] = 0  # Underflow
                        elif i_reco == nx_reco - 1:
                            coords[0] = axis_reco_x.GetNbins() + 1  # Overflow
                        else:
                            coords[0] = i_reco  # Normal bin
                        
                        # Similar for other dimensions...
                        if j_reco == 0:
                            coords[1] = 0
                        elif j_reco == ny_reco - 1:
                            coords[1] = axis_reco_y.GetNbins() + 1
                        else:
                            coords[1] = j_reco
                            
                        if i_gen == 0:
                            coords[2] = 0
                        elif i_gen == nx_gen - 1:
                            coords[2] = axis_gen_x.GetNbins() + 1
                        else:
                            coords[2] = i_gen
                            
                        if j_gen == 0:
                            coords[3] = 0
                        elif j_gen == ny_gen - 1:
                            coords[3] = axis_gen_y.GetNbins() + 1
                        else:
                            coords[3] = j_gen
                    else:
                        coords[0] = i_reco
                        coords[1] = j_reco
                        coords[2] = i_gen
                        coords[3] = j_gen
                    
                    try:
                        content = thnf_4d.GetBinContent(coords)
                        if content > 0:
                            response_2d.SetBinContent(flat_reco_bin, flat_gen_bin, content)
                            filled_bins += 1
                            total_content += content
                            
                            # Track overflow content
                            if include_overflow:
                                is_overflow = (i_reco == 0 or i_reco == nx_reco-1 or 
                                             j_reco == 0 or j_reco == ny_reco-1 or
                                             i_gen == 0 or i_gen == nx_gen-1 or
                                             j_gen == 0 or j_gen == ny_gen-1)
                                if is_overflow:
                                    overflow_content += content
                        else:
                            response_2d.SetBinContent(flat_reco_bin, flat_gen_bin, 0)
                    except:
                        continue
        
        if i_reco % 20 == 0:
            current_max = reco_end - 1 if include_overflow else nx_reco
            print(f"  Processed reco slice {i_reco}/{current_max}")
    
    print(f"Filled {filled_bins} response matrix elements")
    print(f"Total response content: {total_content:.6f}")
    
    if include_overflow:
        print(f"Overflow/underflow content: {overflow_content:.6f} ({100*overflow_content/total_content:.1f}%)")
    
    # Create mapping information
    mapping = {
        'reco_nx': nx_reco, 'reco_ny': ny_reco,
        'gen_nx': nx_gen, 'gen_ny': ny_gen,
        'include_overflow': include_overflow,
        'reco_x_min': axis_reco_x.GetBinLowEdge(1),
        'reco_x_max': axis_reco_x.GetBinUpEdge(axis_reco_x.GetNbins()),
        'reco_y_min': axis_reco_y.GetBinLowEdge(1),
        'reco_y_max': axis_reco_y.GetBinUpEdge(axis_reco_y.GetNbins()),
        'gen_x_min': axis_gen_x.GetBinLowEdge(1),
        'gen_x_max': axis_gen_x.GetBinUpEdge(axis_gen_x.GetNbins()),
        'gen_y_min': axis_gen_y.GetBinLowEdge(1),
        'gen_y_max': axis_gen_y.GetBinUpEdge(axis_gen_y.GetNbins())
    }
    
    return response_2d, mapping

def create_2d_response_matrix_from_thnf(thnf_4d, name_suffix="_response_2d"):
    """
    SAFE version: Convert 4D THnF to 2D TH2 response matrix without segfaults.
    Uses explicit bin iteration instead of global bin enumeration.
    """
    
    print("Converting 4D THnF to 2D response matrix (SAFE VERSION)...")
    
    if thnf_4d.GetNdimensions() != 4:
        raise ValueError(f"Expected 4D THnF, got {thnf_4d.GetNdimensions()}D")
    
    # Get dimensions
    axis_reco_x = thnf_4d.GetAxis(0)
    axis_reco_y = thnf_4d.GetAxis(1)
    axis_gen_x = thnf_4d.GetAxis(2)
    axis_gen_y = thnf_4d.GetAxis(3)
    
    nx_reco = axis_reco_x.GetNbins()
    ny_reco = axis_reco_y.GetNbins()
    nx_gen = axis_gen_x.GetNbins()
    ny_gen = axis_gen_y.GetNbins()
    
    total_reco_bins = nx_reco * ny_reco
    total_gen_bins = nx_gen * ny_gen
    
    print(f"THnF dimensions: reco({nx_reco}x{ny_reco}) -> gen({nx_gen}x{ny_gen})")
    print(f"Flattened response matrix: {total_reco_bins} x {total_gen_bins}")
    
    # Create 2D response matrix
    response_name = thnf_4d.GetName() + name_suffix
    response_2d = ROOT.TH2D(
        response_name, "2D Response Matrix",
        total_reco_bins, 0.5, total_reco_bins + 0.5,  # Reco (X-axis)
        total_gen_bins, 0.5, total_gen_bins + 0.5      # Gen (Y-axis)
    )
    
    # SAFE ITERATION: Use explicit bin loops instead of global bin enumeration
    print("Filling 2D response matrix using SAFE iteration...")
    
    filled_bins = 0
    coords = array('i', [0, 0, 0, 0])
    total_content = 0.0
    
    # Iterate through all possible bin combinations explicitly
    for i_reco in range(1, nx_reco + 1):
        for j_reco in range(1, ny_reco + 1):
            # Reco bin in flattened space (column-major)
            #flat_reco_bin = (j_reco - 1) * nx_reco + i_reco
            flat_reco_bin = (i_reco - 1) * ny_reco + j_reco
            
            for i_gen in range(1, nx_gen + 1):
                for j_gen in range(1, ny_gen + 1):
                    # Gen bin in flattened space (column-major)
                    #flat_gen_bin = (j_gen - 1) * nx_gen + i_gen
                    flat_gen_bin = (i_gen - 1) * ny_gen + j_gen
                    
                    coords[0] = i_reco
                    coords[1] = j_reco
                    coords[2] = i_gen
                    coords[3] = j_gen
                    
                    try:
                        content = thnf_4d.GetBinContent(coords)
                        if content > 0:
                            response_2d.SetBinContent(flat_reco_bin, flat_gen_bin, content)
                            filled_bins += 1
                            total_content += content
                        else:
                            response_2d.SetBinContent(flat_reco_bin, flat_gen_bin, 0)
                    except:
                        continue
        
        if i_reco % 20 == 0:
            print(f"  Processed reco slice {i_reco}/{nx_reco}")
    
    print(f"Filled {filled_bins} response matrix elements")
    print(f"Total response content: {total_content:.6f}")
    
    # Create mapping information
    mapping = {
        'reco_nx': nx_reco, 'reco_ny': ny_reco,
        'gen_nx': nx_gen, 'gen_ny': ny_gen,
        'reco_x_min': axis_reco_x.GetBinLowEdge(1),
        'reco_x_max': axis_reco_x.GetBinUpEdge(nx_reco),
        'reco_y_min': axis_reco_y.GetBinLowEdge(1),
        'reco_y_max': axis_reco_y.GetBinUpEdge(ny_reco),
        'gen_x_min': axis_gen_x.GetBinLowEdge(1),
        'gen_x_max': axis_gen_x.GetBinUpEdge(nx_gen),
        'gen_y_min': axis_gen_y.GetBinLowEdge(1),
        'gen_y_max': axis_gen_y.GetBinUpEdge(ny_gen)
    }
    
    return response_2d, mapping


def handle_roounfold_covariance(unfolder):
    """
    Simple covariance handling - just convert matrix to 2D histogram and keep it flat.
    No unflattening needed!
    """
    
    try:
        cov_obj = unfolder.Eunfold()
        
        if cov_obj is None:
            print("No covariance matrix returned")
            return None, 'none'
        
        # If it's already a 2D histogram, return as-is
        if hasattr(cov_obj, 'GetNbinsX') and hasattr(cov_obj, 'GetNbinsY'):
            print(f"Covariance is 2D histogram: {cov_obj.GetNbinsX()}x{cov_obj.GetNbinsY()}")
            return cov_obj, 'histogram_2d'
        
        # If it's a matrix, convert to 2D histogram
        elif hasattr(cov_obj, 'GetNrows'):
            print(f"Converting matrix to 2D histogram: {cov_obj.GetNrows()}x{cov_obj.GetNcols()}")
            
            nrows = cov_obj.GetNrows()
            ncols = cov_obj.GetNcols()
            
            # Create 2D covariance histogram in flattened bin space
            cov_hist = ROOT.TH2D(
                "covariance_matrix", "Covariance Matrix",
                nrows, 0.5, nrows + 0.5,  # Flattened bin indices
                ncols, 0.5, ncols + 0.5   # Flattened bin indices
            )
            
            # Copy matrix values to histogram
            for i in range(nrows):
                for j in range(ncols):
                    cov_value = cov_obj(i, j)
                    cov_hist.SetBinContent(i + 1, j + 1, cov_value)
            
            print(f"Created covariance histogram: {cov_hist.GetNbinsX()}x{cov_hist.GetNbinsY()}")
            return cov_hist, 'matrix_converted'
        
        # If it's a 1D histogram (diagonal only), convert to 2D diagonal
        elif hasattr(cov_obj, 'GetNbinsX'):
            print(f"Converting 1D diagonal to 2D: {cov_obj.GetNbinsX()} bins")
            
            nbins = cov_obj.GetNbinsX()
            
            # Create diagonal 2D covariance matrix
            cov_hist = ROOT.TH2D(
                "covariance_diagonal", "Diagonal Covariance Matrix",
                nbins, 0.5, nbins + 0.5,
                nbins, 0.5, nbins + 0.5
            )
            
            # Fill diagonal elements only
            for i in range(1, nbins + 1):
                error = cov_obj.GetBinContent(i)
                variance = error * error  # Convert error to variance
                cov_hist.SetBinContent(i, i, variance)
            
            print(f"Created diagonal covariance: {cov_hist.GetNbinsX()}x{cov_hist.GetNbinsY()}")
            return cov_hist, 'diagonal_converted'
        
        else:
            print(f"Unknown covariance type: {type(cov_obj)}")
            return None, 'unknown'
    
    except Exception as e:
        print(f"Error handling covariance: {e}")
        return None, 'failed'

def verify_binning_preservation(original_hist, reconstructed_hist):
    """
    Verify that bin edges are correctly preserved.
    """
    print("\n=== BINNING VERIFICATION ===")
    
    # Check dimensions
    if (original_hist.GetNbinsX() != reconstructed_hist.GetNbinsX() or 
        original_hist.GetNbinsY() != reconstructed_hist.GetNbinsY()):
        print("❌ ERROR: Different number of bins!")
        return False
    
    nx = original_hist.GetNbinsX()
    ny = original_hist.GetNbinsY()
    
    # Check X-axis bin edges
    x_edges_match = True
    for i in range(1, nx + 2):  # Include upper edge of last bin
        orig_edge = original_hist.GetXaxis().GetBinLowEdge(i)
        reco_edge = reconstructed_hist.GetXaxis().GetBinLowEdge(i)
        
        if abs(orig_edge - reco_edge) > 1e-15:
            print(f"❌ X-axis bin {i}: original={orig_edge:.6f}, reconstructed={reco_edge:.6f}")
            x_edges_match = False
            break
    
    # Check Y-axis bin edges  
    y_edges_match = True
    for j in range(1, ny + 2):  # Include upper edge of last bin
        orig_edge = original_hist.GetYaxis().GetBinLowEdge(j)
        reco_edge = reconstructed_hist.GetYaxis().GetBinLowEdge(j)
        
        if abs(orig_edge - reco_edge) > 1e-15:
            print(f"❌ Y-axis bin {j}: original={orig_edge:.6f}, reconstructed={reco_edge:.6f}")
            y_edges_match = False
            break
    
    if x_edges_match and y_edges_match:
        print("✅ All bin edges preserved correctly!")
        
        # Show a few sample bin edges to confirm
        print("Sample X bin edges:")
        for i in [1, 2, 3, nx//2, nx]:
            orig = original_hist.GetXaxis().GetBinLowEdge(i)
            reco = reconstructed_hist.GetXaxis().GetBinLowEdge(i)
            print(f"  Bin {i}: {orig:.6f} (both)")
            
        print("Sample Y bin edges:")
        for j in [1, 2, 3, ny//2, ny]:
            orig = original_hist.GetYaxis().GetBinLowEdge(j)
            reco = reconstructed_hist.GetYaxis().GetBinLowEdge(j)
            print(f"  Bin {j}: {orig:.6f} (both)")
            
        return True
    else:
        print("❌ Bin edges do not match!")
        return False

def extract_diagonal_errors_1d(cov_hist_2d):
    """
    Error extraction
    """
    
    if cov_hist_2d is None:
        return None
    
    nbins = min(cov_hist_2d.GetNbinsX(), cov_hist_2d.GetNbinsY())
    
    print(f"\nDEBUGGING ERROR EXTRACTION:")
    print(f"Covariance matrix: {nbins}x{nbins}")
    
    # Check a few diagonal elements
    print("Sample diagonal elements (should be variances):")
    for i in [1, 10, 50, 100, min(nbins, 500)]:
        if i <= nbins:
            variance = np.float64(cov_hist_2d.GetBinContent(i, i))
            error = np.sqrt(np.maximum(0.0, variance))
            print(f"  Bin {i}: variance={variance:.6f}, sqrt(variance)={error:.6f}")
    
    # Create error histogram
    error_hist = ROOT.TH1D(
        "errors_1d_debug", "Per-bin Errors (Debug)",
        nbins, 0.5, nbins + 0.5
    )
    
    # Check if we're dealing with errors vs variances
    total_diagonal = 0
    for i in range(1, nbins + 1):
        variance = cov_hist_2d.GetBinContent(i, i)
        total_diagonal += variance
        
        # Extract error (square root of variance)
        error = np.sqrt(np.maximum(0.0, variance))
        error_hist.SetBinContent(i, error)
    
    print(f"Total diagonal sum: {total_diagonal:.6f}")
    print(f"Average diagonal: {total_diagonal/nbins:.6f}")
    print(f"Error histogram integral: {error_hist.Integral():.6f}")
    
    # Check if diagonal elements look like errors^2 (variances) or errors
    sample_diag = cov_hist_2d.GetBinContent(10, 10)
    print(f"Sample diagonal element: {sample_diag:.6f}")
    
    if sample_diag < 1e-6:
        print("⚠️  Diagonal elements are very small - might be correct")
    elif sample_diag > 1.0:
        print("⚠️  Diagonal elements are large - check if these are variances")
    
    return error_hist

def extract_nominal_errors_1d(error_1d_full, gen_mapping, name_suffix="_nominal"):
    """
    Extract nominal range errors from the full error histogram that includes overflow bins.
    """
    
    if error_1d_full is None:
        return None
    
    include_overflow = gen_mapping.get('include_overflow', False)
    
    if not include_overflow:
        # No overflow, return as-is
        return error_1d_full
    
    # Get dimensions
    nx_full = gen_mapping['nx']  # Includes overflow
    ny_full = gen_mapping['ny']  # Includes overflow
    nx_nominal = nx_full - 2     # Nominal bins only
    ny_nominal = ny_full - 2     # Nominal bins only
    
    total_nominal_bins = nx_nominal * ny_nominal
    
    print(f"Extracting nominal errors from {error_1d_full.GetNbinsX()} -> {total_nominal_bins} bins")
    
    # Create new error histogram for nominal bins only
    error_name = error_1d_full.GetName() + name_suffix
    error_1d_nominal = ROOT.TH1D(
        error_name, "Nominal Per-bin Errors",
        total_nominal_bins, 0.5, total_nominal_bins + 0.5
    )
    
    # Extract nominal errors using the same mapping as covariance
    for i in range(1, nx_nominal + 1):  # 1-based nominal bins
        for j in range(1, ny_nominal + 1):
            # In the full matrix, nominal bins start at i=1, j=1 (after underflow i=0, j=0)
            i_full = i + 1  # Skip underflow bin
            j_full = j + 1  # Skip underflow bin
            
            # Convert to flat bin indices
            full_flat_bin = (i_full - 1) * ny_full + j_full + 1  # ROOT 1-based
            nominal_flat_bin = (i - 1) * ny_nominal + j           # ROOT 1-based
            
            # Extract error
            error_value = error_1d_full.GetBinContent(full_flat_bin)
            error_1d_nominal.SetBinContent(nominal_flat_bin, error_value)
    
    print(f"Extracted nominal errors: integral = {error_1d_nominal.Integral():.6f}")
    
    return error_1d_nominal
    """
    Debug version of error extraction
    """
    
    if cov_hist_2d is None:
        return None
    
    nbins = min(cov_hist_2d.GetNbinsX(), cov_hist_2d.GetNbinsY())
    
    print(f"\nDEBUGGING ERROR EXTRACTION:")
    print(f"Covariance matrix: {nbins}x{nbins}")
    
    # Check a few diagonal elements
    print("Sample diagonal elements (should be variances):")
    for i in [1, 10, 50, 100, min(nbins, 500)]:
        if i <= nbins:
            variance = np.float64(cov_hist_2d.GetBinContent(i, i))
            error = np.sqrt(np.maximum(0.0, variance))
            print(f"  Bin {i}: variance={variance:.6f}, sqrt(variance)={error:.6f}")
    
    # Create error histogram
    error_hist = ROOT.TH1D(
        "errors_1d_debug", "Per-bin Errors (Debug)",
        nbins, 0.5, nbins + 0.5
    )
    
    # Check if we're dealing with errors vs variances
    total_diagonal = 0
    for i in range(1, nbins + 1):
        variance = cov_hist_2d.GetBinContent(i, i)
        total_diagonal += variance
        
        # Extract error (square root of variance)
        error = np.sqrt(np.maximum(0.0, variance))
        error_hist.SetBinContent(i, error)
    
    print(f"Total diagonal sum: {total_diagonal:.6f}")
    print(f"Average diagonal: {total_diagonal/nbins:.6f}")
    print(f"Error histogram integral: {error_hist.Integral():.6f}")
    
    # Check if diagonal elements look like errors^2 (variances) or errors
    sample_diag = cov_hist_2d.GetBinContent(10, 10)
    print(f"Sample diagonal element: {sample_diag:.6f}")
    
    if sample_diag < 1e-6:
        print("⚠️  Diagonal elements are very small - might be correct")
    elif sample_diag > 1.0:
        print("⚠️  Diagonal elements are large - check if these are variances")
    
    return error_hist

def extract_nominal_covariance_matrix(cov_2d_full, gen_mapping, name_suffix="_nominal"):
    """
    Extract the nominal range sub-matrix from the full covariance matrix that includes overflow bins.
    """
    
    if cov_2d_full is None:
        return None
    
    include_overflow = gen_mapping.get('include_overflow', False)
    
    if not include_overflow:
        # No overflow, return as-is
        print("No overflow bins, returning covariance matrix as-is")
        return cov_2d_full
    
    # Get dimensions
    nx_full = gen_mapping['nx']  # Includes overflow
    ny_full = gen_mapping['ny']  # Includes overflow
    nx_nominal = nx_full - 2     # Nominal bins only
    ny_nominal = ny_full - 2     # Nominal bins only
    
    total_full_bins = nx_full * ny_full
    total_nominal_bins = nx_nominal * ny_nominal
    
    print(f"Extracting nominal covariance matrix:")
    print(f"  Full matrix: {total_full_bins}x{total_full_bins} (includes overflow)")
    print(f"  Nominal matrix: {total_nominal_bins}x{total_nominal_bins} (nominal only)")
    
    # Create new covariance matrix for nominal bins only
    cov_name = cov_2d_full.GetName() + name_suffix
    cov_2d_nominal = ROOT.TH2D(
        cov_name, "Nominal Covariance Matrix",
        total_nominal_bins, 0.5, total_nominal_bins + 0.5,
        total_nominal_bins, 0.5, total_nominal_bins + 0.5
    )
    
    # Map nominal bins from full matrix to nominal matrix
    nominal_bin_map = []  # Maps nominal_flat_bin -> full_flat_bin
    
    for i in range(1, nx_nominal + 1):  # 1-based nominal bins
        for j in range(1, ny_nominal + 1):
            # In the full matrix, nominal bins start at i=1, j=1 (after underflow i=0, j=0)
            i_full = i + 1  # Skip underflow bin
            j_full = j + 1  # Skip underflow bin
            
            # Convert to flat bin indices
            full_flat_bin = (i_full - 1) * ny_full + j_full  # 0-based then +1 for ROOT
            nominal_bin_map.append(full_flat_bin)
    
    print(f"Created mapping for {len(nominal_bin_map)} nominal bins")
    
    # Extract the nominal submatrix
    extracted_elements = 0
    total_variance = 0.0
    
    for nom_i, full_i in enumerate(nominal_bin_map):
        for nom_j, full_j in enumerate(nominal_bin_map):
            # Get covariance element from full matrix
            cov_value = cov_2d_full.GetBinContent(full_i + 1, full_j + 1)  # ROOT is 1-based
            
            # Set in nominal matrix
            cov_2d_nominal.SetBinContent(nom_i + 1, nom_j + 1, cov_value)
            
            if nom_i == nom_j:  # Diagonal element
                total_variance += cov_value
            
            if cov_value != 0:
                extracted_elements += 1
    
    print(f"Extracted {extracted_elements} non-zero covariance elements")
    print(f"Total variance (diagonal sum): {total_variance:.6f}")
    
    # Verify dimensions match the unfolded result
    expected_size = nx_nominal * ny_nominal
    actual_size = cov_2d_nominal.GetNbinsX()
    
    if actual_size == expected_size:
        print(f"✅ Covariance matrix dimensions match unfolded result: {expected_size}x{expected_size}")
    else:
        print(f"❌ ERROR: Dimension mismatch! Expected {expected_size}, got {actual_size}")
    
    return cov_2d_nominal

def check_unfolding_sanity(data_hist, unfolded_hist, name=""):
    """Check if unfolding results make sense."""
    if unfolded_hist is None:
        print(f"ERROR: {name} - Unfolded histogram is None!")
        return False
        
    print(f"\n=== SANITY CHECK: {name} ===")
    
    data_integral = data_hist.Integral()
    unfolded_integral = unfolded_hist.Integral()
    ratio = np.float64(unfolded_integral) / np.float64(data_integral) if data_integral > 0 else 0
    
    print(f"Data integral: {data_integral:.2f}")
    print(f"Unfolded integral: {unfolded_integral:.2f}")
    print(f"Ratio: {ratio:.3f}")
    
    # Check for negative bins
    negative_bins = 0
    total_bins = unfolded_hist.GetNbinsX() * unfolded_hist.GetNbinsY()
    
    for i in range(1, unfolded_hist.GetNbinsX() + 1):
        for j in range(1, unfolded_hist.GetNbinsY() + 1):
            if unfolded_hist.GetBinContent(i, j) < 0:
                negative_bins += 1
    
    print(f"Negative bins: {negative_bins}/{total_bins}")
    
    if 0.8 <= ratio <= 1.2 and negative_bins < total_bins * 0.05:
        print("✅ Unfolding looks GOOD")
        return True
    else:
        print("⚠️  Unfolding looks PROBLEMATIC")
        return False

def check_covariance_matrix_sanity(cov_2d, data_1d, name=""):
    """
    Check if the covariance matrix makes sense.
    """
    print(f"\n=== COVARIANCE MATRIX SANITY CHECK {name} ===")
    
    if cov_2d is None:
        print("No covariance matrix to check")
        return False
    
    # Check dimensions
    nx = cov_2d.GetNbinsX()
    ny = cov_2d.GetNbinsY()
    expected_size = data_1d.GetNbinsX()
    
    print(f"Covariance dimensions: {nx}x{ny}")
    print(f"Expected size: {expected_size}x{expected_size}")
    
    if nx != expected_size or ny != expected_size:
        print(f"❌ Size mismatch! Expected {expected_size}x{expected_size}")
        return False
    
    # Check if matrix is symmetric (it should be)
    asymmetric_elements = 0
    for i in range(1, min(11, nx + 1)):  # Check first 10x10 block
        for j in range(1, min(11, ny + 1)):
            val_ij = cov_2d.GetBinContent(i, j)
            val_ji = cov_2d.GetBinContent(j, i)
            
            if abs(val_ij - val_ji) > 1e-10:
                asymmetric_elements += 1
    
    print(f"Asymmetric elements in 10x10 block: {asymmetric_elements}/100")
    
    # Check diagonal vs off-diagonal
    diagonal_sum = 0.0
    off_diagonal_sum = 0.0
    
    for i in range(1, min(101, nx + 1)):  # Check first 100x100 block
        for j in range(1, min(101, ny + 1)):
            val = cov_2d.GetBinContent(i, j)
            if i == j:
                diagonal_sum += abs(val)
            else:
                off_diagonal_sum += abs(val)
    
    if diagonal_sum > 0:
        off_diag_fraction = off_diagonal_sum / diagonal_sum
        print(f"Off-diagonal/diagonal ratio: {off_diag_fraction:.3f}")
        
        if off_diag_fraction > 10:
            print("⚠️  WARNING: Very large off-diagonal elements")
    
    # Check for negative diagonal elements (shouldn't happen)
    negative_diagonal = 0
    zero_diagonal = 0
    
    for i in range(1, nx + 1):
        val = cov_2d.GetBinContent(i, i)
        if val < 0:
            negative_diagonal += 1
        elif val == 0:
            zero_diagonal += 1
    
    print(f"Negative diagonal elements: {negative_diagonal}/{nx}")
    print(f"Zero diagonal elements: {zero_diagonal}/{nx}")
    
    if negative_diagonal > 0:
        print("❌ ERROR: Negative diagonal elements in covariance matrix!")
        return False
    
    return True

def perform_2d_unfolding_via_flattening_with_overflow(reco_2d, gen_2d, response_thnf, data_2d, 
                                                     n_iterations=4, method="Bayes", 
                                                     include_overflow=True):
    """
    2D EEC unfolding that includes overflow in response matrix but returns only nominal range.
    
    This approach:
    1. Captures all migration patterns (including to/from overflow) in response matrix
    2. Unfolds using full information 
    3. Returns only the nominal range in final result (overflow events are "lost")
    """
    
    print("="*60)
    print(f"PERFORMING 2D UNFOLDING WITH OVERFLOW CAPTURE")
    print(f"Response matrix includes overflow: {include_overflow}")
    print(f"Final result: nominal range only")
    print("="*60)
    
    # Step 1: Flatten with overflow option for response matrix accuracy
    print("\nStep 1: Flattening 2D histograms...")
    reco_1d, reco_mapping = flatten_2d_histogram_with_overflow(reco_2d, "_reco_flat", include_overflow)
    gen_1d, gen_mapping = flatten_2d_histogram_with_overflow(gen_2d, "_gen_flat", include_overflow)
    data_1d, data_mapping = flatten_2d_histogram_with_overflow(data_2d, "_data_flat", include_overflow)
    
    # Step 2: Create response matrix with overflow to capture all migrations
    print(f"\nStep 2: Creating response matrix with overflow capture...")
    response_2d, response_mapping = create_2d_response_matrix_from_thnf_with_overflow(
        response_thnf, include_overflow=include_overflow)
    
    # Step 3: RooUnfold setup
    print("\nStep 3: Creating RooUnfoldResponse...")
    roounfold_response = ROOT.RooUnfoldResponse(
        reco_1d, gen_1d, response_2d,
        "response_2d_overflow", "2D Response with Overflow Capture", False
    )

    # Step 4: Perform unfolding with full overflow-aware response matrix
    if method == "Bayes":
        print(f"\nStep 4: Performing {method} unfolding with {n_iterations} iterations...")
        unfolder = ROOT.RooUnfoldBayes(roounfold_response, data_1d, n_iterations)
    elif method == "BinByBin":
        print(f"\nStep 4: Performing {method} unfolding...")
        unfolder = ROOT.RooUnfoldBinByBin(roounfold_response, data_1d)
    elif method == "Invert":
        print(f"\nStep 4: Performing {method} unfolding...")
        unfolder = ROOT.RooUnfoldInvert(roounfold_response, data_1d)
    else:
        raise ValueError(f"Unknown method: {method}")

    unfolded_1d = unfolder.Hunfold()
    
    print(f"1D unfolding completed")
    print(f"  Data integral: {data_1d.Integral():.6f}")
    print(f"  Unfolded integral: {unfolded_1d.Integral():.6f}")
    
    # Step 5: Enhanced covariance handling with debugging
    print("\nStep 5: Extracting and debugging covariance matrix...")
    cov_2d_full, cov_type = handle_roounfold_covariance(unfolder)
    
    # Debug full covariance matrix
    cov_ok = check_covariance_matrix_sanity(cov_2d_full, data_1d, "(Full with Overflow)")
    
    if cov_2d_full is not None and cov_ok:
        print(f"Full covariance matrix: {cov_2d_full.GetNbinsX()}x{cov_2d_full.GetNbinsY()} (includes overflow)")
        
        # Extract diagonal errors from full matrix
        error_1d_full = extract_diagonal_errors_1d(cov_2d_full)
        
        # Extract ONLY nominal range from both covariance and errors
        print("\nExtracting nominal range from full covariance matrix...")
        cov_2d = extract_nominal_covariance_matrix(cov_2d_full, gen_mapping, "_nominal")
        error_1d_nominal = extract_nominal_errors_1d(error_1d_full, gen_mapping, "_nominal")
        
        # Unflatten the nominal 1D errors to 2D (this should now match perfectly)
        if error_1d_nominal is not None:
            # Create a nominal mapping for unflattening
            nominal_mapping = gen_mapping.copy()
            if include_overflow:
                nominal_mapping['nx'] = gen_mapping['nx'] - 2
                nominal_mapping['ny'] = gen_mapping['ny'] - 2
                nominal_mapping['include_overflow'] = False
                # Keep the nominal edges (already extracted during flattening)
                if 'x_edges' in gen_mapping:
                    nominal_mapping['x_edges'] = gen_mapping['x_edges'][1:-1]  # Skip overflow edges
                if 'y_edges' in gen_mapping:
                    nominal_mapping['y_edges'] = gen_mapping['y_edges'][1:-1]  # Skip overflow edges
            else:
                nominal_mapping = gen_mapping
            
            error_2d = unflatten_1d_to_2d(error_1d_nominal, nominal_mapping, "_error_2d")
            print(f"Error 2D integral after unflattening: {error_2d.Integral():.6f}")
        else:
            error_2d = None
    else:
        print("❌ Covariance matrix issues - no error extraction")
        error_2d = None
        cov_2d = None
    
    # Step 6: Unflatten to nominal range only (ignore overflow bins)
    print("\nStep 6: Unflattening to nominal range (discarding overflow)...")
    unfolded_2d = unflatten_1d_to_2d_nominal_only(unfolded_1d, gen_mapping, "_unfolded_2d")
    
    print(f"\n=== UNFOLDING SUMMARY ===")
    print(f"Input data: {data_2d.Integral():.6f}")
    print(f"Unfolded (nominal): {unfolded_2d.Integral():.6f}")
    print(f"Efficiency (kept in nominal): {unfolded_2d.Integral()/data_2d.Integral():.3f}")
    
    return unfolded_2d, error_2d, cov_2d, response_2d

def perform_2d_unfolding_via_flattening(reco_2d, gen_2d, response_thnf, data_2d, 
                                        n_iterations=4, method="Bayes"):
    """
    Enhanced version with detailed error debugging.
    """
    
    print("="*60)
    print("PERFORMING 2D UNFOLDING WITH ERROR DEBUG")
    print("="*60)
    
    # Steps 1-4: Same as before
    print("\nStep 1: Flattening 2D histograms...")
    reco_1d, reco_mapping = flatten_2d_histogram(reco_2d, "_reco_flat")
    gen_1d, gen_mapping = flatten_2d_histogram(gen_2d, "_gen_flat")
    data_1d, data_mapping = flatten_2d_histogram(data_2d, "_data_flat")
    
    print(f"\nStep 2: Creating response matrix...")
    response_2d, response_mapping = create_2d_response_matrix_from_thnf(response_thnf)
    
    print("\nStep 3: Creating RooUnfoldResponse...")
    roounfold_response = ROOT.RooUnfoldResponse(
        reco_1d, gen_1d, response_2d,
        "response_2d_debug", "2D Response Debug", False
    )

    print(f"\nStep 4: Performing {method} unfolding with {n_iterations} iterations...")
    if method == "Bayes":
        unfolder = ROOT.RooUnfoldBayes(roounfold_response, data_1d, n_iterations)
    elif method == "BinByBin":
        unfolder = ROOT.RooUnfoldBinByBin(roounfold_response, data_1d)
    elif method == "Invert":
        unfolder = ROOT.RooUnfoldInvert(roounfold_response, data_1d)
    else:
        raise ValueError(f"Unknown method: {method}")

    unfolded_1d = unfolder.Hunfold()
    
    print(f"1D unfolding completed")
    print(f"  Data integral: {data_1d.Integral():.6f}")
    print(f"  Unfolded integral: {unfolded_1d.Integral():.6f}")
    print(f"  Ratio: {unfolded_1d.Integral()/data_1d.Integral():.6f}")
    
    # Step 5: Enhanced covariance handling with debugging
    print("\nStep 5: Extracting and debugging covariance matrix...")
    cov_2d, cov_type = handle_roounfold_covariance(unfolder)
    
    # Debug covariance matrix
    cov_ok = check_covariance_matrix_sanity(cov_2d, data_1d, "(Flattened)")
    
    if cov_2d is not None and cov_ok:
        print(f"Covariance matrix: {cov_2d.GetNbinsX()}x{cov_2d.GetNbinsY()} (flattened space)")
        
        # Extract diagonal errors with debugging
        error_1d = extract_diagonal_errors_1d(cov_2d)
        
        # Unflatten the 1D errors to 2D
        if error_1d is not None:
            error_2d = unflatten_1d_to_2d(error_1d, gen_mapping, "_error_2d")
            print(f"Error 2D integral after unflattening: {error_2d.Integral():.6f}")
        else:
            error_2d = None
    else:
        print("❌ Covariance matrix issues - no error extraction")
        error_2d = None
    
    # Step 6: Unflatten result
    print("\nStep 6: Unflattening result back to 2D...")
    unfolded_2d = unflatten_1d_to_2d(unfolded_1d, gen_mapping, "_unfolded_2d")
    
    print(f"2D unfolding completed!")
    print(f"  Final 2D integral: {unfolded_2d.Integral():.6f}")
    
    # Final debugging
    if error_2d is not None:
        print(f"\n=== FINAL ERROR CHECK ===")
        print(f"Unfolded 2D integral: {unfolded_2d.Integral():.6f}")
        print(f"Error 2D integral: {error_2d.Integral():.6f}")
        print(f"Average relative error: {error_2d.Integral()/unfolded_2d.Integral():.6f}")
        
        # Check a few bins for sanity
        print("Sample relative errors:")
        for i, j in [(10, 5), (50, 10), (100, 15)]:
            if i <= unfolded_2d.GetNbinsX() and j <= unfolded_2d.GetNbinsY():
                result = unfolded_2d.GetBinContent(i, j)
                error = error_2d.GetBinContent(i, j)
                rel_err = error / result if result > 0 else 0
                print(f"  Bin ({i},{j}): result={result:.2f}, error={error:.2f}, rel_err={rel_err:.3f}")
    
    return unfolded_2d, error_2d, cov_2d, response_2d

def apply_efficiency_correction_preserving_relative_errors(hUnf2d, eff_corr):
    """
    Apply efficiency correction while preserving the relative error pattern from unfolding
    """
    
    print("Applying efficiency correction while preserving relative errors...")
    
    for i in range(1, hUnf2d.GetNbinsX() + 1):
        for j in range(1, hUnf2d.GetNbinsY() + 1):
            # Get current values
            content = hUnf2d.GetBinContent(i, j)
            error = hUnf2d.GetBinError(i, j)
            eff_factor = eff_corr.GetBinContent(i, j)
            
            # Calculate relative error BEFORE correction
            rel_error = error / content if content > 0 else 0
            
            # Apply efficiency correction to content only
            new_content = content * eff_factor
            
            # Calculate new error to preserve the same relative error
            new_error = new_content * rel_error
            
            # Set the corrected values
            hUnf2d.SetBinContent(i, j, new_content)
            hUnf2d.SetBinError(i, j, new_error)
    
    print("✅ Efficiency correction applied with preserved relative errors")

def unfold_2d_eec_with_overflow(reco_2d_hist, gen_2d_hist, response_4d_thnf, data_2d_hist, 
                                n_iterations=4, unfolding_method="Invert", include_overflow=True):
    """
    2D EEC unfolding that includes overflow in response matrix but returns only nominal range.
    
    Returns:
    --------
    unfolded_2d : ROOT.TH2D
        Unfolded result in original 2D coordinates (nominal range only)
    error_2d : ROOT.TH2D  
        Per-bin errors in original 2D coordinates (optional, can be None)
    cov_2d : ROOT.TH2D
        Covariance matrix in flattened coordinate space
    """
    
    try:
        unfolded_2d, error_2d, cov_2d, resp_2d = perform_2d_unfolding_via_flattening_with_overflow(
            reco_2d_hist, gen_2d_hist, response_4d_thnf, data_2d_hist, 
            n_iterations, unfolding_method, include_overflow
        )
    
        return unfolded_2d, error_2d, cov_2d, resp_2d
        
    except Exception as e:
        print(f"ERROR in unfold_2d_eec_with_overflow: {e}")
        return None, None, None, None

def unfold_2d_eec(reco_2d_hist, gen_2d_hist, response_4d_thnf, data_2d_hist, 
                  n_iterations=4, unfolding_method="Invert"):
    """
    Simple 2D EEC unfolding that keeps covariance in flattened 2D space.
    
    Returns:
    --------
    unfolded_2d : ROOT.TH2D
        Unfolded result in original 2D coordinates
    error_2d : ROOT.TH2D  
        Per-bin errors in original 2D coordinates (optional, can be None)
    cov_2d : ROOT.TH2D
        Covariance matrix in flattened coordinate space
    """
    
    try:
        unfolded_2d, error_2d, cov_2d, resp_2d = perform_2d_unfolding_via_flattening(
            reco_2d_hist, gen_2d_hist, response_4d_thnf, data_2d_hist, 
            n_iterations, unfolding_method
        )
    
        return unfolded_2d, error_2d, cov_2d, resp_2d
        
    except Exception as e:
        print(f"ERROR in unfold_2d_eec: {e}")
        return None, None, None, None

if __name__ == '__main__':


    ## response matrix
    ## v12 best fine binning with double precision 
    ## v11 coarser r and eij bins than v9
    ## v10 coarser r bin than v9
    ## v9 angular new eij and r bin
    ## v8 angular new eij bin
    ## v7 correspondence table
    ## v6 angular

    
    parser = argparse.ArgumentParser(description='Quick TUnfold Test with Fixed Tau')
    parser.add_argument('--method', type=str, default="Bayes",
                       help='Fixed unfolding method to use (default: Bayes)')
    parser.add_argument('--niter', type=int, default=4,
                       help='number of iterations')
    parser.add_argument('--overflow', action='store_true', default=True,
                       help='Include overflow/underflow bins in response matrix')
    parser.add_argument('--jacobian', type=str, default="r",
                       help='Unfold r or z distribution')
    args = parser.parse_args()

    if args.method == "Bayes":
        surfix = f'_niter{args.niter}'
    else:
        surfix = ''
    
    if args.overflow:
        surfix += '_overflow'

    surfix += f'_{args.jacobian}'
    

    #filenamein = 'response_kk2f4146_qqpy_91.25_v14.root'
    #filenamein = 'response_kk2f4146_qqardcy_91.25_v1.root'
    #filenamein = 'response_apacic105_91.25_v1.root'
    #filenamein = 'response_qqps_91.25_v1.root'
    
    #datafile = 'h_kk2f4146_qqpy_91.25_v5.root'

    #filenameout = f'unfolded_qqpy_91.25_qqardcy_v2{surfix}.root'
    #filenameout = f'unfolded_qqpy_91.25_apacic_v2{surfix}.root'

    #data2dname = f'EEC2d_{args.jacobian}'

    # response v15 -> data v3 vvv tight angular
    # response v14 -> data v1 angular (nominal)
    # response v13 -> data v2 correspondence table
    filenamein = 'response_kk2f4146_qqpy_91.25_v15.root'
    #filenamein = 'response_kk2f4146_qqardcy_91.25_v1.root'
    #filenamein = 'response_apacic105_91.25_v1.root'
    #filenamein = 'response_pythia8_91.25_v2.root'
    #filenamein = 'response_pythia8_dire_91.25_v1.root'
    datafile = 'h_94c_v5.root'
    #datafile = 'h_pythia8_v3.root'
    filenameout = f'unfolded_data_kk2f4146_qqpy_91.25_v3{surfix}.root'
    #filenameout = f'unfolded_data_kk2f4146_qqardcy_91.25_v1{surfix}.root'
    #filenameout = f'unfolded_data_pythia8_dire_91.25_v1{surfix}.root'
    #filenameout = f'unfolded_pythia8_91.25_v1{surfix}.root'
    data2dname = f'EEC2d_{args.jacobian}'
    
    fin = ROOT.TFile.Open(filenamein,'r')
    
    counter = fin.Get("counter").Clone("N")
    _response = fin.Get(f"response2d_eij_{args.jacobian}")

    _reco2d = fin.Get(f'reco2d_eij_{args.jacobian}_match')
    _gen2d = fin.Get(f"gen2d_eij_{args.jacobian}_match")

    gen2d = fin.Get(f'gen2d_eij_{args.jacobian}')
    eff_corr = gen2d.Clone("eff_corr")
    eff_corr.Divide(_gen2d)
    
    reco2d = fin.Get(f'reco2d_eij_{args.jacobian}')
    fake_corr = reco2d.Clone("fake_corr")
    fake_corr.Add(_reco2d, -1)
    fake_corr2 = fake_corr.Clone('fake_corr2')
    print("DEBUG1:", fake_corr.Integral())
    fake_corr.Divide(reco2d)
    print("DEBUG2:", fake_corr.Integral())

    fdata = ROOT.TFile.Open(datafile,'r')
    data2d = fdata.Get(data2dname).Clone('data')

    data2d2 = fdata.Get(data2dname).Clone('data')

    normalization = fin.Get("counter").GetBinContent(2)
    try:
        n = fdata.Get('counter').GetBinContent(2)
    except:
        n = fdata.Get('N').GetBinContent(2)

    print("Norms:", normalization, n, normalization/n)
    data2d.Scale(float(normalization)/n)
    print("DEBUG3:", data2d.Integral())
    
    data2d2.Scale(float(normalization)/n)
    data2d2.Add(fake_corr2, -1)
    
    fake_sub=data2d.Clone("fake_sub")
    apply_efficiency_correction_preserving_relative_errors(fake_sub, fake_corr)
    #fake_sub.Multiply(fake_corr)
    print("DEBUG4:", fake_sub.Integral())
    data2d.Add(fake_sub, -1)
    print("DEBUG5:", data2d.Integral())
    print("DEBUG5:", data2d2.Integral())

    # Choose unfolding approach based on overflow flag
    if args.overflow:
        hUnf2d, hErr2d, hCov2d, hResp2d = unfold_2d_eec_with_overflow(
            _reco2d, _gen2d, _response, data2d,
            n_iterations=args.niter,
            unfolding_method=args.method,
            include_overflow=True
        )
    else:
        hUnf2d, hErr2d, hCov2d, hResp2d = unfold_2d_eec(
            _reco2d, _gen2d, _response, data2d,
            n_iterations=args.niter,
            unfolding_method=args.method
        )

    if hUnf2d is None:
        print("FATAL: Unfolding failed!")
        sys.exit(1)

    # Apply unfolding errors BEFORE efficiency correction
    if hErr2d is not None:
        print("Applying unfolding errors to hUnf2d...")
        
        if (hUnf2d.GetNbinsX() == hErr2d.GetNbinsX() and 
            hUnf2d.GetNbinsY() == hErr2d.GetNbinsY()):
            
            for i in range(1, hUnf2d.GetNbinsX() + 1):
                for j in range(1, hUnf2d.GetNbinsY() + 1):
                    error = hErr2d.GetBinContent(i, j)
                    hUnf2d.SetBinError(i, j, error)
            
            print(f"✅ Applied unfolding errors to {hUnf2d.GetNbinsX()}x{hUnf2d.GetNbinsY()} bins")
            
        else:
            print(f"❌ ERROR: Dimension mismatch between hUnf2d and hErr2d!")
            sys.exit(1)
    else:
        print("⚠️  WARNING: No error histogram returned")
    
    check_unfolding_sanity(data2d, hUnf2d, f"Unfolding Method: {args.method} (overflow={args.overflow})")
    
    # Verify binning preservation
    success = verify_binning_preservation(_reco2d, hUnf2d)
    if not success:
        print("WARNING: Binning was not preserved correctly!")
    else:
        print("✅ Non-uniform binning preserved correctly in unfolded result")
    
    reco = []
    gen = []
    unfold = []

    fout = ROOT.TFile(filenameout, 'recreate')
    fout.cd()

    hResp2d.Write("response_2d")
    hUnf2d.Write("unfolded_2d_noeff")

    apply_efficiency_correction_preserving_relative_errors(hUnf2d, eff_corr)
    
    eec_reco = Proj2D_X(data2d, eijbins[1], eijbins[-1], f"RECO_EEC")
    eec_gen = Proj2D_X(gen2d, eijbins[1], eijbins[-1], f"GEN_EEC")
    eec_unfold = Proj2D_X(hUnf2d, eijbins[1], eijbins[-1], f"UNFOLD_EEC")
    eec_reco.SetDirectory(0)
    eec_gen.SetDirectory(0)
    eec_unfold.SetDirectory(0)

    # In your file writing:
    hUnf2d.Write("unfolded_2d")

    if hErr2d is not None:
        hErr2d.Write("errors_2d")      # Per-bin errors in original 2D space
    
    if hCov2d is not None:
        hCov2d.Write("covariance_2d")  # Covariance matrix in flattened space

    counter.Write("N")
        
    gen2d.Write("gen_2d")
    reco2d.Write("reco_2d")
    data2d.Write("data_2d")

    fake_corr.Write("fake_corr")
    eff_corr.Write("eff_corr")

    _gen2d.Write("gen_2d_match")
    _reco2d.Write("reco_2d_match")

    eec_reco.Write()
    eec_gen.Write()
    eec_unfold.Write()

    fout.Close()
    
    print(f"\n=== UNFOLDING COMPLETE ===")
    print(f"Method: {args.method}")
    if args.method == "Bayes":
        print(f"Iterations: {args.niter}")
    print(f"Overflow handling: {args.overflow}")
    print(f"Output file: {filenameout}")
    print(f"Final unfolded integral: {hUnf2d.Integral():.6f}")
