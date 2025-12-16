import React, { useState } from 'react';
import {
  BeakerIcon,
  CubeIcon,
  PlayIcon,
  ClockIcon,
  CheckCircleIcon,
  ExclamationTriangleIcon,
  InformationCircleIcon,
  ArrowDownTrayIcon,
  ChartBarIcon
} from '@heroicons/react/24/outline';
import { clsx } from 'clsx';

const EnhancedPredictions = () => {
  const [inputType, setInputType] = useState('smiles');
  const [inputValue, setInputValue] = useState('');
  const [selectedModels, setSelectedModels] = useState(['clintox']);
  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState(null);

  const availableModels = [
    {
      id: 'tox21',
      name: 'Tox21',
      type: 'Toxicity Assessment',
      description: 'Multi-target toxicity prediction across 12 biological pathways',
      outputType: 'classification',
      accuracy: '89.5%',
      category: 'toxicity',
      examples: 'Thalidomide, Warfarin, Aspirin'
    },
    {
      id: 'clintox',
      name: 'Clinical Toxicity',
      type: 'Toxicity Assessment',
      description: 'FDA approval and clinical trial toxicity risk prediction',
      outputType: 'classification',
      accuracy: '94.2%',
      category: 'toxicity',
      examples: 'Ibuprofen, Acetaminophen, Diclofenac'
    },
    {
      id: 'bbbp',
      name: 'Blood-Brain Barrier Penetration',
      type: 'ADMET Property (Distribution)',
      description: 'BBB permeability for CNS targeting (critical for neurological drugs)',
      outputType: 'classification',
      accuracy: '91.8%',
      category: 'admet',
      examples: 'Caffeine, Morphine, Dopamine'
    },
    {
      id: 'caco2',
      name: 'Caco-2 Permeability',
      type: 'ADMET Property (Absorption)',
      description: 'Intestinal absorption prediction via epithelial cell permeability',
      outputType: 'regression',
      accuracy: 'R²: 0.87',
      category: 'admet',
      examples: 'Metformin, Atenolol, Propranolol'
    },
    {
      id: 'clearance',
      name: 'Intrinsic Clearance',
      type: 'ADMET Property (Metabolism)',
      description: 'Enzyme-mediated clearance rate (affects dosing frequency and half-life)',
      outputType: 'regression',
      accuracy: 'R²: 0.83',
      category: 'admet',
      examples: 'Midazolam, Propofol, Lidocaine'
    },
    {
      id: 'hlm_clint',
      name: 'HLM Intrinsic Clearance',
      type: 'ADMET Property (Metabolism)',
      description: 'Human liver microsomal clearance (hepatic metabolism rate prediction)',
      outputType: 'regression',
      accuracy: 'R²: 0.85',
      category: 'admet',
      examples: 'Diazepam, Verapamil, Terfenadine'
    }
  ];

  const handleModelToggle = (modelId) => {
    setSelectedModels(prev => {
      if (prev.includes(modelId)) {
        return prev.filter(id => id !== modelId);
      } else {
        return [...prev, modelId];
      }
    });
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!inputValue.trim() || selectedModels.length === 0) return;

    setIsLoading(true);
    setResults(null);

    try {
      const response = await fetch('http://localhost:5000/api/predict', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles: inputValue,
          models: selectedModels
        }),
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const data = await response.json();
      
      if (data.error) {
        throw new Error(data.error);
      }

      setResults(data);

    } catch (error) {
      console.error('Prediction error:', error);
      setResults({
        error: error.message,
        smiles: inputValue
      });
    } finally {
      setIsLoading(false);
    }
  };

  const getInterpretation = (model, prediction, probability) => {
    if (model.outputType === 'classification') {
      if (probability > 0.7) return { text: 'High Risk', color: 'red' };
      if (probability > 0.4) return { text: 'Moderate Risk', color: 'yellow' };
      return { text: 'Low Risk', color: 'green' };
    } else {
      // Regression models
      if (prediction > 0.7) return { text: 'High', color: 'green' };
      if (prediction > 0.3) return { text: 'Moderate', color: 'yellow' };
      return { text: 'Low', color: 'red' };
    }
  };

  const exportResults = (format) => {
    if (!results) return;

    if (format === 'csv') {
      const csv = [
        ['Model', 'Prediction', 'Probability/Value', 'Interpretation'],
        ...Object.entries(results.predictions || {}).map(([modelId, pred]) => {
          const model = availableModels.find(m => m.id === modelId);
          const interp = getInterpretation(model, pred.prediction, pred.probability);
          return [
            model?.name || modelId,
            pred.prediction,
            pred.probability || pred.value || 'N/A',
            interp.text
          ];
        })
      ].map(row => row.join(',')).join('\n');

      const blob = new Blob([csv], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'prediction_results.csv';
      a.click();
    } else if (format === 'json') {
      const json = JSON.stringify(results, null, 2);
      const blob = new Blob([json], { type: 'application/json' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'prediction_results.json';
      a.click();
    }
  };

  return (
    <div className="space-y-6">
      {/* Page Header */}
      <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg shadow-primary-500/10">
        <h1 className="text-2xl font-bold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent mb-2">Multi-Property Prediction</h1>
        <p className="text-gray-400">
          Predict individual ADMET properties and toxicity using independently trained GIN models
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Left Column - Input & Model Selection */}
        <div className="lg:col-span-2 space-y-6">
          {/* Input Section */}
          <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg">
            <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent mb-4">Input Method</h2>
            
            <div className="flex space-x-2 mb-6">
              <button
                onClick={() => setInputType('smiles')}
                className={clsx(
                  'px-4 py-2 rounded-lg font-medium transition-all',
                  inputType === 'smiles'
                    ? 'bg-gradient-to-r from-primary-600 to-accent-600 text-white shadow-lg shadow-primary-500/30'
                    : 'bg-gray-800 text-gray-300 hover:bg-gray-700'
                )}
              >
                SMILES Input
              </button>
              <button
                onClick={() => setInputType('batch')}
                className={clsx(
                  'px-4 py-2 rounded-lg font-medium transition-all',
                  inputType === 'batch'
                    ? 'bg-gradient-to-r from-primary-600 to-accent-600 text-white shadow-lg shadow-primary-500/30'
                    : 'bg-gray-800 text-gray-300 hover:bg-gray-700'
                )}
              >
                Batch Upload
              </button>
            </div>

            {inputType === 'smiles' && (
              <div>
                <label className="block text-sm font-medium text-gray-400 mb-2">
                  SMILES String
                </label>
                <input
                  type="text"
                  value={inputValue}
                  onChange={(e) => setInputValue(e.target.value)}
                  placeholder="Enter SMILES notation (e.g., CC(C)Cc1ccc(cc1)C(C)C(O)=O)"
                  className="w-full px-4 py-3 bg-gray-800 border border-gray-700 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500 font-mono text-sm text-gray-200 placeholder-gray-500"
                />
                <p className="text-xs text-gray-500 mt-2">
                  Enter a valid SMILES string for molecular structure
                </p>
              </div>
            )}

            {inputType === 'batch' && (
              <div className="border-2 border-dashed border-primary-500/30 rounded-lg p-8 text-center bg-gray-800/30 hover:border-primary-500/50 transition-all">
                <BeakerIcon className="h-12 w-12 text-primary-400 mx-auto mb-3" />
                <p className="text-gray-400 mb-2">Upload CSV file with SMILES column</p>
                <button className="px-4 py-2 bg-gradient-to-r from-primary-600 to-accent-600 text-white rounded-lg hover:from-primary-500 hover:to-accent-500 shadow-lg shadow-primary-500/30 transition-all">
                  Browse Files
                </button>
              </div>
            )}
          </div>

          {/* Model Selection */}
          <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent">Select Models</h2>
              <span className="text-sm text-gray-400">
                {selectedModels.length} selected
              </span>
            </div>
            
            <div className="space-y-3">
              {availableModels.map((model) => (
                <div
                  key={model.id}
                  onClick={() => handleModelToggle(model.id)}
                  className={clsx(
                    'p-4 rounded-lg border-2 cursor-pointer transition-all group',
                    selectedModels.includes(model.id)
                      ? 'border-primary-500 bg-gradient-to-br from-primary-500/20 to-accent-500/20 shadow-lg shadow-primary-500/20'
                      : 'border-gray-700 bg-gray-800/30 hover:border-primary-500/50'
                  )}
                >
                  <div className="flex items-start justify-between">
                    <div className="flex items-start space-x-3 flex-1">
                      <div className={clsx(
                        'h-10 w-10 rounded-lg flex items-center justify-center flex-shrink-0 transition-transform group-hover:scale-110',
                        selectedModels.includes(model.id) 
                          ? 'bg-gradient-to-br from-primary-500 to-accent-600 shadow-lg shadow-primary-500/50' 
                          : 'bg-gradient-to-br from-primary-500/20 to-accent-500/20'
                      )}>
                        <CubeIcon className={clsx(
                          'h-5 w-5',
                          selectedModels.includes(model.id) ? 'text-white' : 'text-primary-400'
                        )} />
                      </div>
                      <div className="flex-1">
                        <div className="flex items-center space-x-2">
                          <h3 className="font-semibold text-gray-200">{model.name}</h3>
                          <span className="text-xs px-2 py-0.5 bg-gray-700 text-gray-300 rounded-md">
                            {model.type}
                          </span>
                        </div>
                        <p className="text-sm text-gray-400 mt-1">{model.description}</p>
                        <div className="flex items-center space-x-4 mt-2">
                          <span className="text-xs text-primary-400 font-medium">
                            {model.accuracy}
                          </span>
                        </div>
                        {model.examples && (
                          <div className="mt-2 text-xs text-gray-500">
                            Examples: {model.examples}
                          </div>
                        )}
                      </div>
                    </div>
                    <div className={clsx(
                      'h-5 w-5 rounded border-2 flex items-center justify-center flex-shrink-0',
                      selectedModels.includes(model.id)
                        ? 'border-primary-500 bg-primary-500 shadow-sm shadow-primary-500/50'
                        : 'border-gray-600'
                    )}>
                      {selectedModels.includes(model.id) && (
                        <CheckCircleIcon className="h-4 w-4 text-white" />
                      )}
                    </div>
                  </div>
                </div>
              ))}
            </div>

            <button
              onClick={handleSubmit}
              disabled={!inputValue.trim() || selectedModels.length === 0 || isLoading}
              className="w-full mt-6 px-6 py-3 bg-gradient-to-r from-primary-600 to-accent-600 text-white font-semibold rounded-lg hover:from-primary-500 hover:to-accent-500 disabled:from-gray-700 disabled:to-gray-700 disabled:cursor-not-allowed transition-all flex items-center justify-center space-x-2 shadow-lg shadow-primary-500/30"
            >
              {isLoading ? (
                <>
                  <ClockIcon className="h-5 w-5 animate-spin" />
                  <span>Processing...</span>
                </>
              ) : (
                <>
                  <PlayIcon className="h-5 w-5" />
                  <span>Run Prediction</span>
                </>
              )}
            </button>
          </div>
        </div>

        {/* Right Column - Results */}
        <div className="lg:col-span-1">
          <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 sticky top-6 shadow-lg">
            <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent mb-4">Prediction Results</h2>
            
            {!results && !isLoading && (
              <div className="text-center py-12">
                <ChartBarIcon className="h-16 w-16 text-gray-700 mx-auto mb-3" />
                <p className="text-gray-500 text-sm">
                  Results will appear here after prediction
                </p>
              </div>
            )}

            {isLoading && (
              <div className="text-center py-12">
                <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary-500 mx-auto mb-3 shadow-sm shadow-primary-500/50"></div>
                <p className="text-gray-400">Running models...</p>
              </div>
            )}

            {results && results.error && (
              <div className="bg-red-900/20 border border-red-500/50 rounded-lg p-4">
                <div className="flex items-start space-x-2">
                  <ExclamationTriangleIcon className="h-5 w-5 text-red-400 flex-shrink-0 mt-0.5" />
                  <div>
                    <h3 className="font-medium text-red-300">Prediction Error</h3>
                    <p className="text-sm text-red-400 mt-1">{results.error}</p>
                  </div>
                </div>
              </div>
            )}

            {results && results.predictions && (
              <div className="space-y-4">
                {/* Results Table */}
                <div className="border border-gray-700 rounded-lg overflow-hidden">
                  <table className="w-full text-sm">
                    <thead className="bg-gray-800/50">
                      <tr>
                        <th className="px-3 py-2 text-left text-xs font-medium text-gray-400 uppercase">Model</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-gray-400 uppercase">Result</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-gray-800">
                      {Object.entries(results.predictions).map(([modelId, pred]) => {
                        const model = availableModels.find(m => m.id === modelId);
                        const interpretation = getInterpretation(model, pred.prediction, pred.probability);
                        
                        return (
                          <tr key={modelId} className="hover:bg-gray-800/30 transition-colors">
                            <td className="px-3 py-3">
                              <div className="font-medium text-gray-200">{model?.name}</div>
                              <div className="text-xs text-gray-500">{model?.type}</div>
                            </td>
                            <td className="px-3 py-3">
                              <div className={clsx(
                                'inline-flex px-2 py-1 text-xs font-medium rounded-md',
                                interpretation.color === 'red' && 'bg-red-900/30 text-red-300 border border-red-500/30',
                                interpretation.color === 'yellow' && 'bg-yellow-900/30 text-yellow-300 border border-yellow-500/30',
                                interpretation.color === 'green' && 'bg-green-900/30 text-green-300 border border-green-500/30'
                              )}>
                                {interpretation.text}
                              </div>
                              <div className="text-xs text-gray-500 mt-1">
                                {model?.outputType === 'classification' 
                                  ? `Confidence: ${(pred.probability * 100).toFixed(1)}%`
                                  : `Value: ${pred.prediction.toFixed(3)}`
                                }
                              </div>
                            </td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>

                <div className="bg-primary-500/10 border border-primary-500/30 rounded-lg p-3 backdrop-blur-sm">
                  <div className="flex items-start space-x-2">
                    <InformationCircleIcon className="h-5 w-5 text-primary-400 flex-shrink-0 mt-0.5" />
                    <p className="text-xs text-gray-300">
                      Predictions based on trained GIN models. Results should be validated experimentally.
                    </p>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default EnhancedPredictions;
