
from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
import tempfile
import os
import torch
from typing import Optional
import asyncio
from datetime import datetime

# Import your existing model
from livae.agent import agent

# Import API models and utilities
from .models import (
    AnnDataSummary, AgentParameters, TrainingConfig, TrainingState, 
    TrainingMetrics, EmbeddingResult, UploadResponse, TrainingResponse,
    ShapeInfo
)
from .utils import (
    load_anndata_from_path, create_anndata_summary, apply_qc_filtering,
    save_embedding_as_csv
)

# Initialize FastAPI app
app = FastAPI(
    title="Single Cell Deep Learning Agent",
    description="LiVAE training interface",
    version="0.1.0"
)

# Add CORS middleware for Next.js frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],  # Next.js default port
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global state (single session)
class AppState:
    def __init__(self):
        self.adata = None
        self.data_summary: Optional[AnnDataSummary] = None
        self.agent_instance = None
        self.training_state = TrainingState()
        self.temp_files = []  # Track temporary files for cleanup
        
    def clear_data(self):
        """Clear current dataset and agent"""
        self.adata = None
        self.data_summary = None
        self.agent_instance = None
        self.training_state = TrainingState()
        self._cleanup_temp_files()
    
    def _cleanup_temp_files(self):
        """Clean up temporary files"""
        for file_path in self.temp_files:
            try:
                if os.path.exists(file_path):
                    os.unlink(file_path)
            except:
                pass
        self.temp_files.clear()

# Global app state
state = AppState()

@app.post("/upload-data", response_model=UploadResponse)
async def upload_data(file: UploadFile = File(...)):
    """
    Upload AnnData file and extract summary information.
    Replaces any currently loaded dataset.
    """
    # Validate file type
    if not file.filename.endswith('.h5ad'):
        raise HTTPException(
            status_code=400, 
            detail="Only .h5ad files are supported"
        )
    
    try:
        # Clear previous data
        state.clear_data()
        
        # Save uploaded file temporarily
        with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as temp_file:
            content = await file.read()
            temp_file.write(content)
            temp_file_path = temp_file.name
        
        state.temp_files.append(temp_file_path)
        
        # Load AnnData and create summary
        state.adata = load_anndata_from_path(temp_file_path)
        state.data_summary = create_anndata_summary(
            state.adata, 
            filename=file.filename,
            file_size_bytes=len(content)
        )
        
        return UploadResponse(
            success=True,
            message=f"Successfully loaded {file.filename}",
            data_summary=state.data_summary
        )
        
    except Exception as e:
        state.clear_data()
        raise HTTPException(
            status_code=500, 
            detail=f"Failed to process file: {str(e)}"
        )

@app.get("/data-summary", response_model=AnnDataSummary)
async def get_data_summary():
    """Get summary of currently loaded dataset"""
    if state.data_summary is None:
        raise HTTPException(
            status_code=404, 
            detail="No dataset loaded"
        )
    
    return state.data_summary

@app.post("/start-training", response_model=TrainingResponse)
async def start_training(
    parameters: AgentParameters,
    config: TrainingConfig,
    background_tasks: BackgroundTasks
):
    """Start training on the loaded dataset"""
    if state.adata is None:
        raise HTTPException(
            status_code=400, 
            detail="No dataset loaded. Upload data first."
        )
    
    if state.training_state.is_running:
        raise HTTPException(
            status_code=400, 
            detail="Training already in progress"
        )
    
    try:
        # Apply QC filtering if thresholds provided
        training_data = state.adata
        if config.qcparams:
            training_data = apply_qc_filtering(state.adata, config.qcparams)
        
        # Auto-detect device
        device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
        
        # Initialize agent with parameters
        state.agent_instance = agent(
            adata=training_data,
            layer=parameters.layer,
            percent=parameters.percent,
            irecon=parameters.irecon,
            lorentz=parameters.lorentz,
            beta=parameters.beta,
            dip=parameters.dip,
            tc=parameters.tc,
            info=parameters.info,
            hidden_dim=parameters.hidden_dim,
            latent_dim=parameters.latent_dim,
            i_dim=parameters.i_dim,
            lr=parameters.lr,
            device=device
        )
        
        # Initialize training state
        state.training_state = TrainingState(
            is_running=True,
            total_epochs=config.epochs,
            started_at=datetime.now()
        )
        
        # Start training in background
        background_tasks.add_task(run_training, config.epochs)
        
        return TrainingResponse(
            success=True,
            message=f"Training started with {config.epochs} epochs on {device}",
        )
        
    except Exception as e:
        state.training_state = TrainingState()
        raise HTTPException(
            status_code=500, 
            detail=f"Failed to start training: {str(e)}"
        )


async def run_training(epochs: int):
    """Background task for training execution with 10-epoch averaging"""
    try:
        # Run training using your agent's fit method with custom progress tracking
        for epoch in range(epochs):
            if not state.training_state.is_running:
                break  # Allow for early stopping
                
            # Load data and perform training step (mimicking your agent.fit logic)
            data = state.agent_instance.load_data()
            state.agent_instance.step(data)
            
            state.training_state.current_epoch = epoch + 1
            
            # Update metrics every 10 epochs (matching your agent logic)
            if (epoch + 1) % 10 == 0:
                # Calculate averages of the last 10 epochs instead of single values
                avg_metrics = calculate_average_metrics(epoch + 1)
                
                if avg_metrics:  # Only update if we have valid metrics
                    metrics = TrainingMetrics(
                        epoch=epoch + 1,
                        loss=avg_metrics['loss'],
                        ari=avg_metrics['ari'],
                        nmi=avg_metrics['nmi'], 
                        asw=avg_metrics['asw'],
                        ch=avg_metrics['ch'],
                        db=avg_metrics['db'],
                        pc=avg_metrics['pc']
                    )
                    
                    state.training_state.latest_metrics = metrics
                    state.training_state.history.append(metrics)
            
            # Small async sleep to prevent blocking
            await asyncio.sleep(0.001)
        
        # Training completed - calculate final averaged metrics
        final_metrics = calculate_average_metrics(epochs, final=True)
        if final_metrics:
            final_training_metrics = TrainingMetrics(
                epoch=epochs,
                loss=final_metrics['loss'],
                ari=final_metrics['ari'],
                nmi=final_metrics['nmi'], 
                asw=final_metrics['asw'],
                ch=final_metrics['ch'],
                db=final_metrics['db'],
                pc=final_metrics['pc']
            )
            state.training_state.latest_metrics = final_training_metrics
            # Only append if it's different from the last recorded epoch
            if not state.training_state.history or state.training_state.history[-1].epoch != epochs:
                state.training_state.history.append(final_training_metrics)
        
        state.training_state.is_running = False
        state.training_state.completed_at = datetime.now()
        
    except Exception as e:
        state.training_state.is_running = False
        state.training_state.error_message = str(e)


def calculate_average_metrics(current_epoch: int, final: bool = False) -> dict:
    """
    Calculate average metrics over the last 10 epochs
    
    Args:
        current_epoch: Current epoch number
        final: Whether this is the final calculation (use all available data)
        
    Returns:
        Dictionary with averaged metrics or None if insufficient data
    """
    try:
        # Determine how many epochs to average over
        if final:
            # For final calculation, use up to last 10 epochs or all available
            epochs_to_average = min(10, len(state.agent_instance.score))
        else:
            # For regular updates, always use exactly 10 epochs
            epochs_to_average = 10
            
        # Check if we have enough data
        if len(state.agent_instance.score) < epochs_to_average or len(state.agent_instance.loss) < epochs_to_average:
            return None
            
        # Get the last N epochs of scores and losses
        recent_scores = state.agent_instance.score[-epochs_to_average:]
        recent_losses = state.agent_instance.loss[-epochs_to_average:]
        
        # Calculate averages for each metric
        # Scores format: each score is [ari, nmi, asw, ch, db, pc]
        avg_ari = sum(score[0] for score in recent_scores) / len(recent_scores)
        avg_nmi = sum(score[1] for score in recent_scores) / len(recent_scores)
        avg_asw = sum(score[2] for score in recent_scores) / len(recent_scores)
        avg_ch = sum(score[3] for score in recent_scores) / len(recent_scores)
        avg_db = sum(score[4] for score in recent_scores) / len(recent_scores)
        avg_pc = sum(score[5] for score in recent_scores) / len(recent_scores)
        
        # Loss format: each loss is [total_loss, recon_loss, kl_loss, ...]
        avg_loss = sum(loss[0] for loss in recent_losses) / len(recent_losses)
        
        return {
            'loss': float(avg_loss),
            'ari': float(avg_ari),
            'nmi': float(avg_nmi),
            'asw': float(avg_asw),
            'ch': float(avg_ch),
            'db': float(avg_db),
            'pc': float(avg_pc)
        }
        
    except (IndexError, TypeError, ValueError) as e:
        print(f"Error calculating average metrics: {e}")
        return None


@app.get("/training-progress", response_model=TrainingState)
async def get_training_progress():
    """Get current training progress and metrics"""
    return state.training_state

@app.post("/stop-training")
async def stop_training():
    """Stop current training (if running)"""
    if not state.training_state.is_running:
        raise HTTPException(
            status_code=400, 
            detail="No training in progress"
        )
    
    state.training_state.is_running = False
    return {"success": True, "message": "Training stopped"}

@app.get("/embeddings/interpretable", response_model=EmbeddingResult)
async def get_interpretable_embeddings():
    """Get information about interpretable embeddings"""
    if state.agent_instance is None:
        raise HTTPException(
            status_code=400, 
            detail="No trained model available"
        )
    
    try:
        # Extract embedding shape information
        embedding = state.agent_instance.get_iembed()
        shape = ShapeInfo(n_obs=embedding.shape[0], n_vars=embedding.shape[1])
        
        return EmbeddingResult(
            embedding_type="interpretable",
            shape=shape,
            download_ready=True,
            extracted_at=datetime.now()
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500, 
            detail=f"Failed to extract interpretable embeddings: {str(e)}"
        )

@app.get("/embeddings/latent", response_model=EmbeddingResult)
async def get_latent_embeddings():
    """Get information about latent embeddings"""
    if state.agent_instance is None:
        raise HTTPException(
            status_code=400, 
            detail="No trained model available"
        )
    
    try:
        # Extract embedding shape information
        embedding = state.agent_instance.get_latent()
        shape = ShapeInfo(n_obs=embedding.shape[0], n_vars=embedding.shape[1])
        
        return EmbeddingResult(
            embedding_type="latent",
            shape=shape,
            download_ready=True,
            extracted_at=datetime.now()
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500, 
            detail=f"Failed to extract latent embeddings: {str(e)}"
        )

@app.get("/download/embeddings/{embedding_type}")
async def download_embeddings(embedding_type: str):
    """Download embeddings as CSV file"""
    if embedding_type not in ["interpretable", "latent"]:
        raise HTTPException(
            status_code=400,
            detail="Invalid embedding type. Use 'interpretable' or 'latent'"
        )
    
    if state.agent_instance is None:
        raise HTTPException(
            status_code=400, 
            detail="No trained model available"
        )
    
    try:
        # Extract embeddings
        if embedding_type == "interpretable":
            embedding = state.agent_instance.get_iembed()
        else:  # latent
            embedding = state.agent_instance.get_latent()
        
        # Save as CSV
        csv_path = save_embedding_as_csv(embedding, embedding_type)
        
        # Return file download
        return FileResponse(
            path=csv_path,
            filename=os.path.basename(csv_path),
            media_type='text/csv'
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to generate download: {str(e)}"
        )

@app.get("/status")
async def get_app_status():
    """Get overall application status"""
    return {
        "data_loaded": state.adata is not None,
        "model_trained": state.agent_instance is not None,
        "training_running": state.training_state.is_running,
        "device": str(torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu'))
    }


@app.get("/")
async def root():
    """Welcome endpoint"""
    return {
        "message": "Single Cell Deep Learning Agent API", 
        "version": "0.1.0",
        "docs": "http://127.0.0.1:8000/docs",
        "status": "http://127.0.0.1:8000/status"
    }

# Cleanup on shutdown
@app.on_event("shutdown")
async def shutdown_event():
    state.clear_data()
