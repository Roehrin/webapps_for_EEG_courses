/* ==========================================
   MULTI-QUIZ SYSTEM - CORE FUNCTIONS
   Add this to ALL quiz pages
   ========================================== */

// Configuration: Add all your quiz page URLs here
const QUIZ_CONFIG = {
  quizPages: [
    'topography_quiz_multi.html',
    'graph_quiz_multi.html',
    'source_localization_quiz_multi.html',
	'vector_quiz_multi.html'
    // Add more quiz pages as needed
  ],
  
  // Optional: Weight certain quizzes to appear more frequently
  // Default is equal probability for all
  weights: {
    'topography_quiz_multi.html': 2,
    'graph_quiz_multi.html': 2,
    'source_localization_quiz_multi.html': 1,
	'vector_quiz_multi.html':2
  }
};

/* ==========================================
   SCORE MANAGEMENT (Universal)
   ========================================== */

function getScore() {
  const correct = parseInt(sessionStorage.getItem('quizCorrect') || '0');
  const total = parseInt(sessionStorage.getItem('quizTotal') || '0');
  return { correct, total };
}

function updateScoreDisplay() {
  const { correct, total } = getScore();
  const scoreEl = document.getElementById('scoreDisplay');
  if (scoreEl) {
    scoreEl.textContent = `Score: ${correct}/${total}`;
  }
}

function incrementScore(isCorrect) {
  const { correct, total } = getScore();
  sessionStorage.setItem('quizCorrect', isCorrect ? correct + 1 : correct);
  sessionStorage.setItem('quizTotal', total + 1);
  updateScoreDisplay();
}

function resetScore() {
  if (confirm('Are you sure you want to reset your score?')) {
    sessionStorage.removeItem('quizCorrect');
    sessionStorage.removeItem('quizTotal');
    sessionStorage.removeItem('quizHistory');
    updateScoreDisplay();
  }
}

/* ==========================================
   NAVIGATION SYSTEM
   ========================================== */

// Get a random quiz page (truly random - can repeat)
function getRandomQuizPage() {
  const availablePages = QUIZ_CONFIG.quizPages;
  
  // Apply weights if configured
  const weightedPages = [];
  availablePages.forEach(page => {
    const weight = QUIZ_CONFIG.weights[page] || 1;
    for (let i = 0; i < weight; i++) {
      weightedPages.push(page);
    }
  });
  
  const randomIndex = Math.floor(Math.random() * weightedPages.length);
  return weightedPages[randomIndex];
}

// Get current page name
function getCurrentPageName() {
  return window.location.pathname.split('/').pop() || 'index.html';
}

// Navigate to next random quiz
function goToNextQuiz() {
  const nextPage = getRandomQuizPage();
  window.location.href = nextPage;
}

// Track quiz history (optional - prevents immediate repeats)
function addToHistory(pageName) {
  const history = JSON.parse(sessionStorage.getItem('quizHistory') || '[]');
  history.push({
    page: pageName,
    timestamp: Date.now()
  });
  
  // Keep only last 10 entries
  if (history.length > 10) {
    history.shift();
  }
  
  sessionStorage.setItem('quizHistory', JSON.stringify(history));
}

function getQuizHistory() {
  return JSON.parse(sessionStorage.getItem('quizHistory') || '[]');
}

// Smart next quiz (avoids recent repeats if possible)
function goToNextQuizSmart() {
  const history = getQuizHistory();
  const recentPages = history.slice(-3).map(h => h.page);
  
  let availablePages = QUIZ_CONFIG.quizPages.filter(
    page => !recentPages.includes(page)
  );
  
  // If all pages were recent, just use all pages
  if (availablePages.length === 0) {
    availablePages = QUIZ_CONFIG.quizPages;
  }
  
  const randomIndex = Math.floor(Math.random() * availablePages.length);
  const nextPage = availablePages[randomIndex];
  
  addToHistory(getCurrentPageName());
  window.location.href = nextPage;
}

/* ==========================================
   INITIALIZATION
   Call this on page load for each quiz
   ========================================== */

function initMultiQuizSystem() {
  // Update score display
  updateScoreDisplay();
  
  // Make functions globally available
  window.resetScore = resetScore;
  window.goToNextQuiz = goToNextQuiz;
  window.incrementScore = incrementScore;
  window.updateScoreDisplay = updateScoreDisplay;
}

/* ==========================================
   HELPER: Update "Next Question" button
   Call this to modify existing next buttons
   ========================================== */

function setupNextButton(buttonId = 'nextBtn') {
  const btn = document.getElementById(buttonId);
  if (btn) {
    // Remove existing onclick
    btn.removeAttribute('onclick');
    
    // Add new navigation
    btn.addEventListener('click', goToNextQuiz);
    
    // Optional: Update button text to be more generic
    btn.textContent = 'Next Question →';
  }
}

// Auto-initialize on DOM load
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', initMultiQuizSystem);
} else {
  initMultiQuizSystem();
}