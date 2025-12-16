import React, { useState, useRef, useEffect } from 'react';
import { 
  PaperAirplaneIcon, 
  UserIcon, 
  BeakerIcon,
  ClipboardDocumentIcon,
  TrashIcon,
  SparklesIcon
} from '@heroicons/react/24/outline';

const Chat = () => {
  const [messages, setMessages] = useState([
    {
      id: 1,
      role: 'assistant',
      content: 'Welcome to the AI Scientific Assistant for Computational Chemistry.\n\nI provide professional, research-grade explanations for:\n\n• Individual ADMET property predictions (Absorption, Distribution, Metabolism, Excretion, Toxicity)\n• Clinical toxicity risk interpretation\n• Pharmacokinetic property analysis\n• Structure-property relationships\n• Each property is predicted by an independently trained GIN model\n\nNote: ADMET is not a single model—each property has its own dedicated predictor.\n\nHow may I assist you today?',
      timestamp: new Date()
    }
  ]);
  const [inputMessage, setInputMessage] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const messagesEndRef = useRef(null);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const exampleQuestions = [
    "What does high BBB penetration mean for CNS drugs?",
    "Explain Caco-2 permeability and oral absorption",
    "What is HLM intrinsic clearance and why does it matter?",
    "How should I interpret clinical toxicity predictions?",
    "What's the difference between intrinsic clearance and HLM clearance?",
    "Why are ADMET properties predicted separately?"
  ];

  const handleSendMessage = async () => {
    if (!inputMessage.trim() || isLoading) return;

    const userMessage = {
      id: Date.now(),
      role: 'user',
      content: inputMessage.trim(),
      timestamp: new Date()
    };

    setMessages(prev => [...prev, userMessage]);
    setInputMessage('');
    setIsLoading(true);

    try {
      const response = await fetch('http://localhost:5000/api/chat/ask', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          message: userMessage.content,
          context: 'scientific_admet'
        })
      });

      if (response.ok) {
        const data = await response.json();
        const assistantMessage = {
          id: Date.now() + 1,
          role: 'assistant',
          content: data.response || 'I apologize, but I encountered an issue processing your request. Please try again.',
          timestamp: new Date()
        };
        setMessages(prev => [...prev, assistantMessage]);
      } else {
        throw new Error('Failed to get response from AI assistant');
      }
    } catch (error) {
      console.error('Chat error:', error);
      const errorMessage = {
        id: Date.now() + 1,
        role: 'assistant',
        content: 'I am currently unable to connect to the AI service. Please ensure the backend server is running and try again.',
        timestamp: new Date()
      };
      setMessages(prev => [...prev, errorMessage]);
    } finally {
      setIsLoading(false);
    }
  };

  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSendMessage();
    }
  };

  const clearChat = () => {
    setMessages([
      {
        id: 1,
        role: 'assistant',
        content: 'Chat cleared. How can I assist you with your ADMET analysis today?',
        timestamp: new Date()
      }
    ]);
  };

  const copyToClipboard = (text) => {
    navigator.clipboard.writeText(text);
  };

  const handleExampleQuestion = (question) => {
    setInputMessage(question);
  };

  return (
    <div className="h-[calc(100vh-8rem)] flex flex-col">
      {/* Header */}
      <div className="bg-gradient-to-br from-gray-900 to-black border-b border-primary-500/20 px-6 py-4 flex items-center justify-between rounded-t-xl shadow-lg">
        <div className="flex items-center space-x-3">
          <div className="h-10 w-10 rounded-lg bg-gradient-to-br from-primary-500 to-accent-600 flex items-center justify-center shadow-lg shadow-primary-500/50">
            <SparklesIcon className="h-5 w-5 text-white" />
          </div>
          <div>
            <h1 className="text-lg font-semibold bg-gradient-to-r from-primary-400 to-accent-400 bg-clip-text text-transparent">AI Research Assistant</h1>
            <p className="text-sm text-gray-400">Scientific explanations and model interpretation</p>
          </div>
        </div>
        <button
          onClick={clearChat}
          className="p-2 text-gray-400 hover:text-gray-200 hover:bg-gray-800 rounded-lg transition-colors"
          title="Clear conversation"
        >
          <TrashIcon className="h-5 w-5" />
        </button>
      </div>

      {/* Messages Container */}
      <div className="flex-1 overflow-y-auto bg-black p-6 space-y-4">
        {messages.map((message) => (
          <div
            key={message.id}
            className={`flex ${message.role === 'user' ? 'justify-end' : 'justify-start'}`}
          >
            <div
              className={`flex items-start space-x-3 max-w-3xl ${
                message.role === 'user' ? 'flex-row-reverse space-x-reverse' : ''
              }`}
            >
              <div
                className={`h-8 w-8 rounded-lg flex items-center justify-center flex-shrink-0 ${
                  message.role === 'user'
                    ? 'bg-gradient-to-br from-primary-500 to-accent-600 shadow-lg shadow-primary-500/50'
                    : 'bg-gradient-to-br from-gray-800 to-gray-700'
                }`}
              >
                {message.role === 'user' ? (
                  <UserIcon className="h-5 w-5 text-white" />
                ) : (
                  <BeakerIcon className="h-5 w-5 text-primary-400" />
                )}
              </div>
              
              <div
                className={`rounded-lg px-4 py-3 ${
                  message.role === 'user'
                    ? 'bg-gradient-to-r from-primary-600 to-accent-600 text-white shadow-lg shadow-primary-500/30'
                    : 'bg-gradient-to-br from-gray-900 to-gray-800 text-gray-200 border border-primary-500/20'
                }`}
              >
                <p className="text-sm whitespace-pre-wrap leading-relaxed">
                  {message.content}
                </p>
                <div className="flex items-center justify-between mt-2">
                  <span
                    className={`text-xs ${
                      message.role === 'user' ? 'text-primary-100' : 'text-gray-500'
                    }`}
                  >
                    {message.timestamp.toLocaleTimeString([], {
                      hour: '2-digit',
                      minute: '2-digit'
                    })}
                  </span>
                  {message.role === 'assistant' && (
                    <button
                      onClick={() => copyToClipboard(message.content)}
                      className="ml-2 text-xs text-gray-400 hover:text-gray-200 transition-colors"
                    >
                      <ClipboardDocumentIcon className="h-4 w-4" />
                    </button>
                  )}
                </div>
              </div>
            </div>
          </div>
        ))}
        
        {isLoading && (
          <div className="flex justify-start">
            <div className="flex items-start space-x-3 max-w-3xl">
              <div className="h-8 w-8 rounded-lg bg-gradient-to-br from-gray-800 to-gray-700 flex items-center justify-center">
                <BeakerIcon className="h-5 w-5 text-primary-400" />
              </div>
              <div className="bg-gradient-to-br from-gray-900 to-gray-800 border border-primary-500/20 rounded-lg px-4 py-3">
                <div className="flex space-x-2">
                  <div className="h-2 w-2 bg-primary-500 rounded-full animate-bounce shadow-sm shadow-primary-500/50"></div>
                  <div className="h-2 w-2 bg-primary-500 rounded-full animate-bounce shadow-sm shadow-primary-500/50" style={{ animationDelay: '0.2s' }}></div>
                  <div className="h-2 w-2 bg-primary-500 rounded-full animate-bounce shadow-sm shadow-primary-500/50" style={{ animationDelay: '0.4s' }}></div>
                </div>
              </div>
            </div>
          </div>
        )}
        
        <div ref={messagesEndRef} />
      </div>

      {/* Example Questions */}
      {messages.length === 1 && (
        <div className="bg-gradient-to-br from-gray-900 to-black border-t border-primary-500/20 px-6 py-4">
          <p className="text-sm font-medium text-gray-400 mb-3">Example questions:</p>
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-2">
            {exampleQuestions.map((question, index) => (
              <button
                key={index}
                onClick={() => handleExampleQuestion(question)}
                className="text-left px-3 py-2 text-sm bg-gray-800/50 hover:bg-gray-700/50 text-gray-300 rounded-lg transition-all border border-gray-700 hover:border-primary-500/50"
              >
                {question}
              </button>
            ))}
          </div>
        </div>
      )}

      {/* Input Area */}
      <div className="bg-gradient-to-br from-gray-900 to-black border-t border-primary-500/20 px-6 py-4 rounded-b-xl shadow-lg">
        <div className="flex items-end space-x-3">
          <textarea
            value={inputMessage}
            onChange={(e) => setInputMessage(e.target.value)}
            onKeyPress={handleKeyPress}
            placeholder="Ask about model predictions, ADMET properties, or drug discovery..."
            className="flex-1 px-4 py-3 bg-gray-800 border border-gray-700 text-gray-200 rounded-lg focus:ring-2 focus:ring-primary-500 focus:border-primary-500 resize-none text-sm placeholder-gray-500"
            rows="2"
            disabled={isLoading}
          />
          <button
            onClick={handleSendMessage}
            disabled={!inputMessage.trim() || isLoading}
            className="px-6 py-3 bg-gradient-to-r from-primary-600 to-accent-600 text-white rounded-lg hover:from-primary-500 hover:to-accent-500 disabled:from-gray-700 disabled:to-gray-700 disabled:cursor-not-allowed transition-all flex items-center space-x-2 shadow-lg shadow-primary-500/30"
          >
            <PaperAirplaneIcon className="h-5 w-5" />
            <span className="hidden sm:inline">Send</span>
          </button>
        </div>
        <p className="text-xs text-gray-500 mt-2">
          This AI assistant provides scientific explanations based on drug discovery knowledge. 
          Responses should be validated with experimental data.
        </p>
      </div>
    </div>
  );
};

export default Chat;
