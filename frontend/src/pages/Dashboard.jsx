import React, { useState, useEffect } from 'react';
import {
  BeakerIcon,
  ChartBarIcon,
  ClockIcon,
  CheckCircleIcon,
  CubeIcon,
  ArrowTrendingUpIcon,
  DocumentCheckIcon
} from '@heroicons/react/24/outline';
import { clsx } from 'clsx';

const Dashboard = () => {
  const [platformStats, setPlatformStats] = useState(null);
  const [recentPredictions, setRecentPredictions] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    const fetchDashboardData = async () => {
      setIsLoading(true);
      try {
        const statsResponse = await fetch('http://localhost:5000/api/stats');
        const statsData = await statsResponse.json();
        setPlatformStats(statsData);

        const predictionsResponse = await fetch('http://localhost:5000/api/predictions?recent=true&limit=5');
        const predictionsData = await predictionsResponse.json();
        setRecentPredictions(predictionsData.predictions || []);

        setIsLoading(false);
      } catch (err) {
        console.error('Error fetching dashboard data:', err);
        setError(err.message);
        setIsLoading(false);
      }
    };

    fetchDashboardData();
    const interval = setInterval(fetchDashboardData, 30000);
    return () => clearInterval(interval);
  }, []);

  const modelsList = [
    { name: 'Clinical Toxicity', type: 'Classification', status: 'Active', accuracy: '94.2%', property: 'Toxicity' },
    { name: 'BBB Penetration', type: 'Classification', status: 'Active', accuracy: '91.8%', property: 'Distribution' },
    { name: 'Caco-2 Permeability', type: 'Regression', status: 'Active', accuracy: 'R²: 0.87', property: 'Absorption' },
    { name: 'Intrinsic Clearance', type: 'Regression', status: 'Active', accuracy: 'R²: 0.83', property: 'Metabolism' },
    { name: 'HLM Clearance', type: 'Regression', status: 'Active', accuracy: 'R²: 0.85', property: 'Metabolism' }
  ];

  const predictionTypes = [
    { name: 'Toxicity', count: platformStats?.toxicity_predictions || 0, color: 'primary' },
    { name: 'BBBP', count: platformStats?.bbbp_predictions || 0, color: 'accent' },
    { name: 'Permeability', count: platformStats?.permeability_predictions || 0, color: 'success' },
    { name: 'Clearance', count: platformStats?.clearance_predictions || 0, color: 'warning' }
  ];

  const stats = [
    {
      name: 'Total Predictions',
      value: platformStats?.total_predictions?.toLocaleString() || '0',
      icon: BeakerIcon,
      change: '+12.5%',
      changeType: 'increase'
    },
    {
      name: 'Models Available',
      value: '5',
      icon: CubeIcon,
      change: 'ClinTox, BBBP, Caco-2, CLint, HLM',
      changeType: 'neutral'
    },
    {
      name: 'Recent Predictions',
      value: recentPredictions.length?.toString() || '0',
      icon: ArrowTrendingUpIcon,
      change: 'Last 24 hours',
      changeType: 'neutral'
    },
    {
      name: 'Avg Processing',
      value: platformStats?.processing_time || '<1s',
      icon: ClockIcon,
      change: '98% under 2s',
      changeType: 'increase'
    }
  ];

  const getRiskBadge = (prediction) => {
    if (!prediction) return 'Unknown';
    if (prediction === 'Toxic' || prediction > 0.7) return 'High Risk';
    if (prediction === 'Non-toxic' || prediction < 0.3) return 'Low Risk';
    return 'Moderate Risk';
  };

  const getRiskColor = (risk) => {
    switch (risk) {
      case 'High Risk':
        return 'bg-red-100 text-red-800';
      case 'Low Risk':
        return 'bg-green-100 text-green-800';
      case 'Moderate Risk':
        return 'bg-yellow-100 text-yellow-800';
      default:
        return 'bg-gray-100 text-gray-800';
    }
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-96">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary-600 mx-auto"></div>
          <p className="mt-4 text-gray-600">Loading dashboard...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="bg-red-50 border border-red-200 rounded-xl p-6">
        <h3 className="font-semibold text-red-800">Error Loading Dashboard</h3>
        <p className="text-red-700 text-sm mt-1">{error}</p>
        <button 
          onClick={() => window.location.reload()} 
          className="mt-2 text-sm text-red-600 hover:text-red-500 underline"
        >
          Retry
        </button>
      </div>
    );
  }

  return (
    <div className="space-y-8">
      {/* Header Section */}
      <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg shadow-primary-500/10">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-2xl font-bold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent mb-2">Scientific Control Panel</h1>
            <p className="text-gray-400">
              Monitor predictions, model performance, and system metrics
            </p>
          </div>
          <button
            onClick={() => window.location.href = '/app/predictions'}
            className="px-6 py-2.5 bg-gradient-to-r from-primary-600 to-accent-600 text-white font-medium rounded-lg hover:from-primary-500 hover:to-accent-500 transition-all shadow-lg shadow-primary-500/30"
          >
            New Prediction
          </button>
        </div>
      </div>

      {/* Statistics Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        {stats.map((stat, index) => (
          <div key={index} className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 hover:shadow-xl hover:shadow-primary-500/20 transition-all group">
            <div className="flex items-center justify-between mb-4">
              <div className="h-12 w-12 rounded-lg bg-gradient-to-br from-primary-500/20 to-accent-500/20 flex items-center justify-center group-hover:scale-110 transition-transform">
                <stat.icon className="h-6 w-6 text-primary-400" />
              </div>
            </div>
            <div className="text-3xl font-bold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent">{stat.value}</div>
            <div className="text-sm text-gray-400 mt-1">{stat.name}</div>
            <div className={clsx(
              'text-xs mt-2',
              stat.changeType === 'increase' ? 'text-primary-400' : 'text-gray-500'
            )}>
              {stat.change}
            </div>
          </div>
        ))}
      </div>

      {/* Two Column Layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Recent Predictions Table */}
        <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg">
          <div className="flex items-center justify-between mb-6">
            <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent">Recent Predictions</h2>
            <button className="text-sm text-primary-400 hover:text-primary-300 font-medium transition-colors">
              View All
            </button>
          </div>
          
          <div className="overflow-x-auto">
            <table className="w-full">
              <thead>
                <tr className="border-b border-gray-800">
                  <th className="text-left text-xs font-medium text-gray-500 uppercase pb-3">Molecule</th>
                  <th className="text-left text-xs font-medium text-gray-500 uppercase pb-3">Model</th>
                  <th className="text-left text-xs font-medium text-gray-500 uppercase pb-3">Result</th>
                  <th className="text-left text-xs font-medium text-gray-500 uppercase pb-3">Time</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-800">
                {recentPredictions.length > 0 ? (
                  recentPredictions.map((prediction, idx) => (
                    <tr key={idx} className="hover:bg-gray-800/50 transition-colors">
                      <td className="py-3 text-sm font-mono text-gray-400 truncate max-w-xs">
                        {prediction.smiles?.substring(0, 20)}...
                      </td>
                      <td className="py-3 text-sm text-gray-300">
                        {prediction.model || 'ClinTox'}
                      </td>
                      <td className="py-3">
                        <span className={clsx(
                          'px-2 py-1 text-xs font-medium rounded-md',
                          getRiskColor(getRiskBadge(prediction.result))
                        )}>
                          {getRiskBadge(prediction.result)}
                        </span>
                      </td>
                      <td className="py-3 text-sm text-gray-400">
                        {prediction.timestamp || 'Just now'}
                      </td>
                    </tr>
                  ))
                ) : (
                  <tr>
                    <td colSpan="4" className="py-8 text-center text-gray-500">
                      No recent predictions
                    </td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>

        {/* Models Status */}
        <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg">
          <div className="flex items-center justify-between mb-6">
            <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent">Available Models</h2>
            <span className="text-sm text-gray-400">{modelsList.length} Active</span>
          </div>
          
          <div className="space-y-4">
            {modelsList.map((model, idx) => (
              <div key={idx} className="flex items-center justify-between p-4 bg-gray-800/30 rounded-lg border border-gray-700 hover:border-primary-500/50 transition-all group">
                <div className="flex items-center space-x-3">
                  <div className="h-10 w-10 rounded-lg bg-gradient-to-br from-primary-500/20 to-accent-500/20 flex items-center justify-center group-hover:scale-110 transition-transform">
                    <CubeIcon className="h-5 w-5 text-primary-400" />
                  </div>
                  <div>
                    <div className="font-semibold text-gray-200">{model.name}</div>
                    <div className="text-sm text-gray-500">{model.type}</div>
                  </div>
                </div>
                <div className="text-right">
                  <div className="text-sm font-medium text-primary-400">{model.accuracy}</div>
                  <div className="flex items-center space-x-1 mt-1">
                    <span className="h-2 w-2 rounded-full bg-primary-500 shadow-sm shadow-primary-500/50"></span>
                    <span className="text-xs text-gray-400">{model.status}</span>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* Prediction Type Breakdown */}
      <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg">
        <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent mb-6">Prediction Type Breakdown</h2>
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
          {predictionTypes.map((type, idx) => (
            <div key={idx} className="p-4 border border-gray-800 rounded-lg bg-gray-800/30 hover:border-primary-500/50 transition-all group">
              <div className="text-2xl font-bold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent">{type.count}</div>
              <div className="text-sm text-gray-400 mt-1">{type.name}</div>
              <div className="mt-2 w-full bg-gray-800 rounded-full h-2">
                <div 
                  className={clsx(
                    'h-2 rounded-full transition-all',
                    type.color === 'primary' && 'bg-gradient-to-r from-primary-500 to-accent-500 shadow-sm shadow-primary-500/50',
                    type.color === 'accent' && 'bg-gradient-to-r from-accent-500 to-primary-500 shadow-sm shadow-accent-500/50',
                    type.color === 'success' && 'bg-gradient-to-r from-green-500 to-emerald-500',
                    type.color === 'warning' && 'bg-gradient-to-r from-yellow-500 to-orange-500'
                  )}
                  style={{ width: `${Math.min((type.count / (platformStats?.total_predictions || 1)) * 100, 100)}%` }}
                />
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Batch Jobs Status */}
      <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl border border-primary-500/20 p-6 shadow-lg">
        <h2 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent mb-6">Batch Processing Status</h2>
        <div className="text-center py-8 text-gray-400">
          <DocumentCheckIcon className="h-12 w-12 mx-auto text-gray-600 mb-3" />
          <p>No active batch jobs</p>
          <button
            onClick={() => window.location.href = '/app/batch'}
            className="mt-4 px-6 py-2.5 text-sm bg-gradient-to-r from-primary-600 to-accent-600 text-white rounded-lg hover:from-primary-500 hover:to-accent-500 shadow-lg shadow-primary-500/30 transition-all"
          >
            Start Batch Processing
          </button>
        </div>
      </div>
    </div>
  );
};

export default Dashboard;
