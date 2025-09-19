// define the submodules

// Re-exporting the CellData struct
pub use crate::cellids::cell_data::CellData;

// Re-exporting the MainClass struct
pub use main_class::CellIds10x;

// Import the submodules
mod cell_data;
mod main_class;