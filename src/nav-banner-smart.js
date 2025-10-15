/**
 * MODULAR NAVIGATION BANNER WITH TABBED DYNAMIC MENU
 * Self-contained navigation component for EEG quiz pages
 * Automatically reads menu structure from index.html with h2 tabs and h3 sections
 * 
 * Usage: Add this script tag to any page:
 * <script src="../src/nav-banner.js"></script>
 * 
 * The banner will automatically inject itself at the top of the page
 */

(function() {
  'use strict';

  // Configuration
  const CONFIG = {
    indexPath: '../index.html',
    brandColor: '#96004B',
    logoUrl: 'https://www.unige.ch/cdn/themes/unige2016/img/unige-logo.svg',
    githubUrl: 'https://github.com/Roehrin/webapps_for_EEG_courses/tree/site',
    copyright: '© 2025 Dr Nicolas Roehri'
  };

  // Create and inject styles
  function injectStyles() {
    const styleId = 'nav-banner-styles';
    if (document.getElementById(styleId)) return;

    const styles = `
      body { 
        margin: 0 !important; 
        padding: 0 !important;
      }
      .nav-banner-spacer { 
        display: block;
        width: 100%;
      }
      .nav-banner {
        position: fixed; top: 0; left: 0; right: 0;
        background: ${CONFIG.brandColor}; color: #f1f1f1; z-index: 99999;
        padding: 15px 20px 20px 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.3);
        font-family: sans-serif;
      }
      .nav-banner-spacer { height: auto; min-height: 100px; }
      .nav-banner-content {
        max-width: 1400px; margin: 0 auto; display: flex; align-items: center;
        justify-content: space-between; flex-wrap: wrap; gap: 10px;
      }
      .nav-banner-left { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; }
      .nav-banner-home-btn, .nav-banner-menu-btn {
        background: rgba(255,255,255,0.2); color: white;
        border: 2px solid rgba(255,255,255,0.4);
        padding: 8px 16px; border-radius: 6px; text-decoration: none;
        font-weight: bold; font-size: 14px; transition: all 0.3s ease;
        display: inline-flex; align-items: center; gap: 6px; cursor: pointer;
      }
      .nav-banner-home-btn:hover, .nav-banner-menu-btn:hover {
        background: rgba(255,255,255,0.3); border-color: rgba(255,255,255,0.6);
        transform: translateY(-2px);
      }
      .nav-banner-menu-btn.active { background: rgba(255,255,255,0.4); }
      .nav-banner-menu-btn:disabled { opacity: 0.5; cursor: not-allowed; }
      .nav-banner-title { font-size: 16px; margin: 0; color: rgba(255,255,255,0.9); }
      .nav-banner-right { display: flex; align-items: center; gap: 15px; }
      .nav-banner-logo { height: 40px; filter: brightness(0) invert(1); opacity: 0.9; }
      .nav-banner-copyright { font-size: 14px; color: rgba(255,255,255,0.9); margin: 0; }
      .nav-banner-dropdown {
        position: absolute; top: 100%; left: 0; right: 0; background: #fff;
        box-shadow: 0 4px 12px rgba(0,0,0,0.2);
        max-height: 0; overflow: hidden; transition: max-height 0.4s ease;
      }
      .nav-banner-dropdown.open { max-height: 600px; overflow-y: auto; }
      .nav-banner-tabs {
        display: flex; gap: 0; background: #f0f0f0; border-bottom: 2px solid ${CONFIG.brandColor};
        padding: 0 20px; overflow-x: auto; scrollbar-width: none; -ms-overflow-style: none;
      }
      .nav-banner-tabs::-webkit-scrollbar { display: none; }
      .nav-banner-tab {
        background: transparent; border: none; color: #555; padding: 12px 24px;
        font-size: 15px; font-weight: 600; cursor: pointer;
        border-bottom: 3px solid transparent; transition: all 0.3s ease;
        white-space: nowrap; position: relative; top: 2px;
      }
      .nav-banner-tab:hover { background: rgba(150, 0, 75, 0.1); color: ${CONFIG.brandColor}; }
      .nav-banner-tab.active {
        color: ${CONFIG.brandColor}; border-bottom-color: ${CONFIG.brandColor}; background: #fff;
      }
      .nav-banner-tab-content { display: none; }
      .nav-banner-tab-content.active { display: block; }
      .nav-banner-dropdown-content {
        max-width: 1400px; margin: 0 auto; padding: 20px;
        display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px;
      }
      .nav-banner-section {
        background: #f9f9f9; border-radius: 8px; padding: 15px;
        border-left: 4px solid ${CONFIG.brandColor};
      }
      .nav-banner-section h3 { margin: 0 0 12px 0; color: ${CONFIG.brandColor}; font-size: 16px; }
      .nav-banner-section a {
        color: #333; text-decoration: none; font-size: 14px; display: block;
        padding: 6px 10px; border-radius: 4px; transition: all 0.2s ease;
      }
      .nav-banner-section a:hover {
        background: ${CONFIG.brandColor}; color: white; transform: translateX(4px);
      }
      /* 🌟 Highlight for current page link */
      .nav-banner-section a.current-page {
        font-weight: bold;
        color: ${CONFIG.brandColor};
        background: rgba(150, 0, 75, 0.08);
      }
      .nav-banner-loading { text-align: center; padding: 20px; color: #666; font-style: italic; }
    `;
    const styleElement = document.createElement('style');
    styleElement.id = styleId;
    styleElement.textContent = styles;
    document.head.appendChild(styleElement);
  }

  // Parse index.html to extract hierarchical menu structure (h2 -> h3 -> links)
  async function parseIndexMenu() {
    try {
      const response = await fetch(CONFIG.indexPath);
      if (!response.ok) throw new Error('Failed to fetch index.html');
      const html = await response.text();
      const parser = new DOMParser();
      const doc = parser.parseFromString(html, 'text/html');
      const menuStructure = {};
      const h2Headers = doc.querySelectorAll('h2');

      h2Headers.forEach(h2 => {
        const tabName = h2.textContent.trim();
        const tabSections = {};
        let currentElement = h2.nextElementSibling;
        let currentH3 = null;

        while (currentElement && currentElement.tagName !== 'H2') {
          if (currentElement.tagName === 'H3') {
            currentH3 = currentElement.textContent.trim();
            tabSections[currentH3] = [];
          } else if (currentH3) {
            const anchors = currentElement.querySelectorAll('a');
            anchors.forEach(a => {
              const href = a.getAttribute('href');
              const text = a.textContent.trim();
              const href2 = "." + href;
              if (href && !href2.startsWith('http') && text) {
                tabSections[currentH3].push({ name: text, url: href2 });
              }
            });
          }
          currentElement = currentElement.nextElementSibling;
        }

        const validSections = Object.entries(tabSections).filter(([_, links]) => links.length > 0);
        if (validSections.length > 0) {
          menuStructure[tabName] = Object.fromEntries(validSections);
        }
      });

      return menuStructure;
    } catch (error) {
      console.error('Error parsing index.html:', error);
      return null;
    }
  }

  // Create banner HTML
  function createBanner() {
    const banner = document.createElement('div');
    banner.className = 'nav-banner';
    banner.innerHTML = `
      <div class="nav-banner-content">
        <div class="nav-banner-left">
          <a href="${CONFIG.indexPath}" class="nav-banner-home-btn">
            <span>←</span><span>Index</span>
          </a>
          <button class="nav-banner-menu-btn" id="nav-menu-toggle" disabled>
            <span>☰</span><span>Loading...</span>
          </button>
          <p class="nav-banner-title">Advanced EEG pages collection</p>
        </div>
        <div class="nav-banner-right">
          <img src="${CONFIG.logoUrl}" alt="University of Geneva" class="nav-banner-logo">
          <p class="nav-banner-copyright">${CONFIG.copyright}</p>
        </div>
      </div>
      <div class="nav-banner-dropdown" id="nav-dropdown">
        <div class="nav-banner-loading">Loading menu...</div>
      </div>
    `;
    return banner;
  }

  // Populate dropdown with tabs and sections
  function populateMenu(menuStructure, activeTabIndex = 0) {
    const dropdown = document.getElementById('nav-dropdown');
    const menuButton = document.getElementById('nav-menu-toggle');
    const currentPage = window.location.pathname.split('/').pop();
    const normalizedCurrentPage = normalizeMultiQuizPage(currentPage);

    if (!menuStructure || Object.keys(menuStructure).length === 0) {
      dropdown.innerHTML = '<div class="nav-banner-loading">No pages found in index.html</div>';
      return;
    }

    const tabsHTML = `
      <div class="nav-banner-tabs">
        ${Object.keys(menuStructure).map((tabName, index) => `
          <button class="nav-banner-tab ${index === activeTabIndex ? 'active' : ''}" data-tab="${index}">
            ${tabName}
          </button>
        `).join('')}
      </div>
    `;

    const tabContentsHTML = Object.entries(menuStructure).map(([tabName, sections], tabIndex) => {
      const sectionsHTML = Object.entries(sections).map(([sectionName, links]) => `
        <div class="nav-banner-section">
          <h3>${sectionName}</h3>
          <ul>
            ${links.map(link => {
              const linkPage = link.url.split('/').pop();
              const normalizedLink = normalizeMultiQuizPage(linkPage);
              const isCurrent = linkPage === currentPage || 
                               normalizedLink === normalizedCurrentPage || 
                               linkPage === normalizedCurrentPage;
              return `
                <li>
                  <a href="${link.url}" class="${isCurrent ? 'current-page' : ''}">
                    ${link.name}
                  </a>
                </li>
              `;
            }).join('')}
          </ul>
        </div>
      `).join('');
      return `
        <div class="nav-banner-tab-content ${tabIndex === activeTabIndex ? 'active' : ''}" data-tab-content="${tabIndex}">
          <div class="nav-banner-dropdown-content">${sectionsHTML}</div>
        </div>
      `;
    }).join('');

    dropdown.innerHTML = tabsHTML + tabContentsHTML;

    if (menuButton) {
      menuButton.disabled = false;
      menuButton.querySelector('span:last-child').textContent = 'Browse Pages';
    }

    setupTabSwitching();
  }

  // Handle tab switching
  function setupTabSwitching() {
    const tabs = document.querySelectorAll('.nav-banner-tab');
    const tabContents = document.querySelectorAll('.nav-banner-tab-content');
    tabs.forEach(tab => {
      tab.addEventListener('click', () => {
        const tabIndex = tab.getAttribute('data-tab');
        tabs.forEach(t => t.classList.remove('active'));
        tabContents.forEach(tc => tc.classList.remove('active'));
        tab.classList.add('active');
        const activeContent = document.querySelector(`[data-tab-content="${tabIndex}"]`);
        if (activeContent) activeContent.classList.add('active');
      });
    });
  }

  // Toggle dropdown menu
  function setupMenuToggle() {
    const toggle = document.getElementById('nav-menu-toggle');
    const dropdown = document.getElementById('nav-dropdown');
    if (toggle && dropdown) {
      toggle.addEventListener('click', () => {
        if (!toggle.disabled) {
          dropdown.classList.toggle('open');
          toggle.classList.toggle('active');
        }
      });
      document.addEventListener('click', (e) => {
        if (!e.target.closest('.nav-banner')) {
          dropdown.classList.remove('open');
          toggle.classList.remove('active');
        }
      });
    }
  }

  // Normalize multi-quiz page names to match the index entry
  function normalizeMultiQuizPage(filename) {
    // Map all three quiz variants to the base name in index.html
    const multiQuizVariants = [
      'topography_quiz_multi.html',
      'graph_quiz_multi.html',
      'source_localization_quiz_multi.html'
    ];
    
    if (multiQuizVariants.includes(filename)) {
      return 'topography_quiz_multi.html'; // Return the canonical one listed in index
    }
    
    return filename;
  }

  // Update spacer height to match banner
  function updateSpacerHeight() {
    const banner = document.querySelector('.nav-banner');
    const spacer = document.querySelector('.nav-banner-spacer');
    if (banner && spacer) {
      const height = banner.offsetHeight;
      spacer.style.height = height + 'px';
      console.log('Banner height updated:', height + 'px');
    }
  }

  // Set up ResizeObserver to automatically track banner size changes
  function setupBannerObserver() {
    const banner = document.querySelector('.nav-banner');
    if (!banner) return;

    const resizeObserver = new ResizeObserver(entries => {
      for (let entry of entries) {
        const computedHeight = banner.offsetHeight;
        const spacer = document.querySelector('.nav-banner-spacer');
        if (spacer) {
          spacer.style.height = computedHeight + 'px';
          console.log('Banner height auto-updated:', computedHeight + 'px');
        }
      }
    });

    resizeObserver.observe(banner);
    
    const content = banner.querySelector('.nav-banner-content');
    if (content) {
      resizeObserver.observe(content);
    }
  }

  // Initialize banner
  async function init() {
    const currentPage = window.location.pathname.split('/').pop();
    if (currentPage === 'index.html' || currentPage === '') return;

    injectStyles();

    if (document.body) {
      const banner = createBanner();
      document.body.insertBefore(banner, document.body.firstChild);
      
      const spacer = document.createElement('div');
      spacer.className = 'nav-banner-spacer';
      document.body.insertBefore(spacer, banner.nextSibling);
      
      setupMenuToggle();

      const menuStructure = await parseIndexMenu();

      // Detect which tab this page belongs to
      let activeTabIndex = 0;
      if (menuStructure) {
        const tabs = Object.entries(menuStructure);
        
        // Normalize current page for multi-quiz variants
        const normalizedPage = normalizeMultiQuizPage(currentPage);
        
        for (let i = 0; i < tabs.length; i++) {
          const [tabName, sections] = tabs[i];
          for (const [section, links] of Object.entries(sections)) {
            // Check if any link matches the current or normalized page
            if (links.some(link => {
              const linkPage = link.url.split('/').pop();
              const normalizedLink = normalizeMultiQuizPage(linkPage);
              return linkPage === currentPage || 
                     normalizedLink === normalizedPage ||
                     linkPage === normalizedPage;
            })) {
              activeTabIndex = i;
              break;
            }
          }
        }
      }

      populateMenu(menuStructure, activeTabIndex);
      
      // Set up automatic height tracking with ResizeObserver
      setupBannerObserver();
      
      // Also update on window resize and load as backup
      window.addEventListener('resize', updateSpacerHeight);
      window.addEventListener('load', updateSpacerHeight);
      
      // Wait for logo image to load before final measurement
      const logo = document.querySelector('.nav-banner-logo');
      if (logo) {
        if (logo.complete) {
          setTimeout(updateSpacerHeight, 100);
        } else {
          logo.addEventListener('load', () => {
            setTimeout(updateSpacerHeight, 100);
          });
        }
      }
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();